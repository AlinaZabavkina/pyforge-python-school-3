FROM python:3
WORKDIR /python
COPY . /python
RUN pip install --no-cache-dir -r requirements.txt
CMD [ "python", "./main.py"]


FROM continuumio/miniconda3
WORKDIR /rdkit
COPY . /rdkit
RUN conda install conda-forge::rdkit


FROM python:3
WORKDIR /fastapi
COPY . /fastapi
RUN pip install --no-cache-dir --upgrade -r requirements.txt
EXPOSE 80
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "80"]