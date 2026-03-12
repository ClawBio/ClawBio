FROM python:3.12-slim
WORKDIR /app
COPY requirements.txt ./requirements.txt
COPY bot/requirements.txt ./bot-requirements.txt
RUN pip install --no-cache-dir -r requirements.txt -r bot-requirements.txt
COPY . .
CMD ["python3", "bot/roboterri.py"]
