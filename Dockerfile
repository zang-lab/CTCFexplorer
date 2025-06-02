FROM python:3.12-slim
WORKDIR /app
ADD static .
ADD templates .
COPY app.py .
COPY requirments.txt .
RUN pip install -r requirments.txt
EXPOSE 5000
ENTRYPOINT ["python3", "app.py"]