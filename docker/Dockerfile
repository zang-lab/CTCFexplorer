FROM python:3.12

RUN apt-get update || true && \
    apt-get install -y procps  && \
    rm -rf /var/lib/apt/lists/*

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

COPY ../requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY ../app.py .
COPY ../templates templates/
COPY ../static static/
COPY --chmod=0755 ../load_data load_data/
COPY ../pipeline pipeline/

EXPOSE 5000

CMD ["python", "app.py"]