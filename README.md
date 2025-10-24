# MGIonetrimmer

This repository stores a script to convert FASTQ headers from MGI format to One Lambda format.

## Como usar

```bash
usage: mgionetrimmer.py [-h] --input INPUT --output OUTPUT

Executa cutadapt e corrige headers dos arquivos FASTQ.

options:
  -h, --help       show this help message and exit
  --input INPUT    Caminho de entrada para os arquivos.
  --output OUTPUT  Diretório onde os arquivos de saída serão armazenados.
```

## opcional: limitar paralelismo

```
python processa_cutadapt_headers.py \
  --input "C:\dados\entrada" \
  --output "C:\dados\saida" \
  --max-workers 2
```
## Observações

Se preferir, você pode setar --max-workers 1 e aumentar per_job_threads (alterando a regra no código) para concentrar tudo em um único cutadapt grande.

Esse script foi especialmente desenvolvido para rodar nos equipamentos DNBSEQ-G99.
