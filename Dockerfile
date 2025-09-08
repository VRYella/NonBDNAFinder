FROM python:3.12-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt ./
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy application code
COPY . ./

# Create non-root user for security
RUN useradd -m -u 1000 nbduser && \
    chown -R nbduser:nbduser /app
USER nbduser

# Expose ports for both Streamlit and API
EXPOSE 8501 8000

# Health checks for both services
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl --fail http://localhost:8501/_stcore/health || exit 1

# Default to Streamlit, but allow API mode via environment variable
ENV SERVICE_MODE=streamlit

# Start script to choose service
COPY <<EOF /app/start.sh
#!/bin/bash
if [ "\$SERVICE_MODE" = "api" ]; then
    echo "🚀 Starting NBDFinder REST API..."
    python api.py
elif [ "\$SERVICE_MODE" = "both" ]; then
    echo "🌟 Starting both Streamlit and API services..."
    python api.py &
    streamlit run app.py --server.port=8501 --server.address=0.0.0.0
else
    echo "🎨 Starting NBDFinder Streamlit App..."
    streamlit run app.py --server.port=8501 --server.address=0.0.0.0
fi
EOF

RUN chmod +x /app/start.sh

ENTRYPOINT ["/app/start.sh"]