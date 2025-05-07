import argparse
import subprocess
import os
import gzip
import shutil
import sys
import logging
import multiprocessing
import re

def setup_logging(output_path):
    log_file = os.path.join(output_path, "processamento.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )

def check_cutadapt_exists():
    if shutil.which("cutadapt") is None:
        logging.error("O executável 'cutadapt' não foi encontrado no PATH.")
        sys.exit(1)

def find_paired_fastqs(input_path):
    files = os.listdir(input_path)
    r1_files = sorted([f for f in files if re.search(r'_R1.*\.fastq\.gz$', f)])
    pairs = []

    for r1 in r1_files:
        r2 = r1.replace('_R1', '_R2')
        if r2 in files:
            pairs.append((os.path.join(input_path, r1), os.path.join(input_path, r2)))
        else:
            logging.warning(f"Arquivo pareado para {r1} não encontrado.")
    
    return pairs

def run_cutadapt(pairs, output_path):
    cpus = str(multiprocessing.cpu_count())
    for r1, r2 in pairs:
        sample_name = os.path.basename(r1).split('_R1')[0]
        out_r1 = os.path.join(output_path, f"{sample_name}_R1_trimmed.fastq.gz")
        out_r2 = os.path.join(output_path, f"{sample_name}_R2_trimmed.fastq.gz")

        cmd = [
            "cutadapt",
            "-b", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "-B", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            "-j", cpus,
            "-o", out_r1,
            "-p", out_r2,
            r1, r2
        ]

        logging.info(f"Executando cutadapt para: {os.path.basename(r1)} / {os.path.basename(r2)}")
        try:
            subprocess.run(cmd, check=True)
            logging.info(f"cutadapt finalizado para {sample_name}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Erro ao executar cutadapt para {sample_name}: {e}")
            sys.exit(1)

def fix_fastq_headers(output_path):
    for filename in os.listdir(output_path):
        file_path = os.path.join(output_path, filename)

        if filename.endswith((".fastq.gz",)):
            logging.info(f"Corrigindo headers em: {file_path} (gzip)")
            tmp_file = file_path + ".tmp"
            with gzip.open(file_path, 'rt') as infile, gzip.open(tmp_file, 'wt') as outfile:
                for i, line in enumerate(infile):
                    if i % 4 == 0 and line.startswith('@'):
                        outfile.write(line.replace('/', ' '))
                    else:
                        outfile.write(line)
            os.replace(tmp_file, file_path)

def main():
    parser = argparse.ArgumentParser(description="Executa cutadapt e corrige headers dos arquivos FASTQ.")
    parser.add_argument("--input", required=True, help="Caminho de entrada para os arquivos.")
    parser.add_argument("--output", required=True, help="Diretório onde os arquivos de saída serão armazenados.")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    setup_logging(args.output)

    logging.info("Iniciando o processamento com cutadapt...")
    check_cutadapt_exists()

    pairs = find_paired_fastqs(args.input)
    if not pairs:
        logging.error("Nenhum par de arquivos FASTQ encontrados.")
        sys.exit(1)

    run_cutadapt(pairs, args.output)
    fix_fastq_headers(args.output)

    logging.info("Processamento concluído com sucesso.")

if __name__ == "__main__":
    main()

