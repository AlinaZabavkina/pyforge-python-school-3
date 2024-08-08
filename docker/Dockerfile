FROM python:3
WORKDIR /python
COPY ../src /python
RUN pip install --no-cache-dir -r requirements.txt
CMD [ "python", "./main.py"]


FROM continuumio/miniconda3
WORKDIR /rdkit
COPY ../src /rdkit
RUN conda install conda-forge::rdkit


FROM python:3
WORKDIR /fastapi
COPY ../src /fastapi
RUN pip install --no-cache-dir --upgrade -r requirements.txt
EXPOSE 80
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "80"]