import argparse
import subprocess
import os
import gzip
import shutil
import sys
import logging
import multiprocessing
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

def setup_logging(output_path: str) -> None:
    os.makedirs(output_path, exist_ok=True)
    log_file = os.path.join(output_path, "processamento.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file, mode="w"), logging.StreamHandler(sys.stdout)],
    )

def check_cutadapt_exists() -> None:
    if shutil.which("cutadapt") is None:
        logging.error("O executável 'cutadapt' não foi encontrado no PATH.")
        sys.exit(1)

def find_paired_fastqs(input_path: str):
    files = os.listdir(input_path)
    # Apenas .fastq.gz com _R1/_R2
    r1_files = sorted([f for f in files if re.search(r"_R1.*\.fastq\.gz$", f)])
    pairs = []
    for r1 in r1_files:
        r2 = r1.replace("_R1", "_R2")
        if r2 in files:
            pairs.append((os.path.join(input_path, r1), os.path.join(input_path, r2)))
        else:
            logging.warning(f"Par não encontrado para: {r1}")
    return pairs

def fix_headers_for_file(gz_path: str) -> None:
    """Descompacta para .fastq temporário, altera apenas headers e recompata para o mesmo .gz."""
    temp_fastq = gz_path[:-3] if gz_path.endswith(".gz") else gz_path + ".tmp.fastq"
    logging.info(f"    [headers] Ajustando: {os.path.basename(gz_path)}")
    try:
        with gzip.open(gz_path, "rt") as fin, open(temp_fastq, "w") as fout:
            for i, line in enumerate(fin):
                if i % 4 == 0 and line.startswith("@"):
                    fout.write(line.replace("/", " "))
                else:
                    fout.write(line)
        # Recompacta no mesmo nome original
        with open(temp_fastq, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    finally:
        # Remove temporário se existir
        try:
            if os.path.exists(temp_fastq):
                os.remove(temp_fastq)
        except Exception as e:
            logging.warning(f"    [headers] Falha ao remover temporário {temp_fastq}: {e}")

def run_cutadapt_for_pair(r1: str, r2: str, output_path: str, per_job_threads: int) -> tuple:
    """
    Executa cutadapt para um par e corrige headers dos dois arquivos de saída.
    Retorna (sample_label, success_bool, error_message_or_None)
    """
    sample_label = os.path.basename(r1).split("_R1")[0]
    out_r1 = os.path.join(output_path, os.path.basename(r1))  # mantém o nome original
    out_r2 = os.path.join(output_path, os.path.basename(r2))  # mantém o nome original

    cmd = [
        "cutadapt",
        "-b", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        "-B", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        "-j", str(max(1, per_job_threads)),
        "-o", out_r1,
        "-p", out_r2,
        r1, r2,
    ]

    logging.info(f"[cutadapt] {sample_label}: iniciando ({os.path.basename(r1)} / {os.path.basename(r2)}) com -j {per_job_threads}")
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"[cutadapt] {sample_label}: finalizado")
    except subprocess.CalledProcessError as e:
        return (sample_label, False, f"Erro cutadapt: {e}")
    except Exception as e:
        return (sample_label, False, f"Erro inesperado: {e}")

    # Corrige headers imediatamente após o trim
    try:
        fix_headers_for_file(out_r1)
        fix_headers_for_file(out_r2)
        logging.info(f"[headers] {sample_label}: concluído")
    except Exception as e:
        return (sample_label, False, f"Erro ajustando headers: {e}")

    return (sample_label, True, None)

def main():
    parser = argparse.ArgumentParser(
        description="Executa cutadapt em pares R1/R2 em paralelo e corrige apenas os headers (mantendo nomes .fastq.gz)."
    )
    parser.add_argument("--input", required=True, help="Diretório de entrada com os FASTQ.gz.")
    parser.add_argument("--output", required=True, help="Diretório de saída para os resultados.")
    parser.add_argument("--max-workers", type=int, default=None,
                        help="Número máximo de pares processados em paralelo (padrão: min(#pares, #CPUs)).")
    args = parser.parse_args()

    setup_logging(args.output)
    logging.info("Iniciando processamento…")
    check_cutadapt_exists()

    pairs = find_paired_fastqs(args.input)
    if not pairs:
        logging.error("Nenhum par R1/R2 .fastq.gz encontrado em --input.")
        sys.exit(1)

    total_cpus = multiprocessing.cpu_count()
    max_workers = args.max_workers if args.max_workers and args.max_workers > 0 else min(len(pairs), total_cpus)

    # Divide os CPUs entre os jobs paralelos para evitar oversubscription
    per_job_threads = max(1, total_cpus // max_workers)

    logging.info(f"Pares encontrados: {len(pairs)} | CPUs: {total_cpus} | workers: {max_workers} | threads/job: {per_job_threads}")

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_cutadapt_for_pair, r1, r2, args.output, per_job_threads) for r1, r2 in pairs]
        for fut in as_completed(futures):
            results.append(fut.result())

    failures = [r for r in results if not r[1]]
    if failures:
        logging.error("Alguns pares falharam:")
        for sample_label, _, msg in failures:
            logging.error(f"  - {sample_label}: {msg}")
        sys.exit(1)

    logging.info("Processamento concluído com sucesso para todos os pares.")

if __name__ == "__main__":
    main()
