FROM jvansan/uvicorn-gunicorn:3.7

LABEL Name=nmr_api Version=0.1

RUN mkdir /usr/share/man/man1/
RUN apt-get update && apt-get install -y openjdk-11-jdk-headless

ENV PYTHONVENV=app
RUN conda install -n $PYTHONVENV -c rdkit rdkit
RUN conda install -n $PYTHONVENV -c conda-forge fastapi

COPY ./app /app/app