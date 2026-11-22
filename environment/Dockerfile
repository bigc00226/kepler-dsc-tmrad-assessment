FROM python:3.13-slim-bookworm
ENV OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
RUN pip install --no-cache-dir numpy==2.5.3 scipy==1.18.1
COPY data /app/data
RUN mkdir -p /app/output
WORKDIR /app
