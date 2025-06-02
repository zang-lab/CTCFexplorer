FROM python:3.12-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Set working directory
WORKDIR /app

# Install dependencies early to cache layers
COPY requirments.txt .
RUN pip install --no-cache-dir -r requirments.txt

# Copy source code
COPY app.py .
COPY templates/ templates/
COPY static/ static/

# Expose Flask default port
EXPOSE 5000

# Use a non-root user if needed (optional hardening for production)
# RUN adduser --disabled-password --gecos '' appuser && chown -R appuser /app
# USER appuser

# Set entrypoint
#CMD ["python3", "app.py"]

ENV FLASK_APP=app.py
CMD ["flask", "run", "--host=0.0.0.0", "--port=5000"]
