# Bloco FFT da estação — um container por rádio.
#
# numpy e pyzmq vêm de wheel manylinux no CPython 3.11: imagem slim, sem
# compilador.

FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

COPY pyproject.toml README.md ./
COPY src/ ./src/

RUN pip install --no-cache-dir -e .

EXPOSE 5582

CMD ["python", "-m", "grs_fft.main", "--help"]
