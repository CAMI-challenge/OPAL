FROM python:3.7.12-slim

ADD *.py /usr/local/bin/
ADD src /usr/local/bin/src
ADD src/utils /usr/local/bin/src/utils
ADD requirements /requirements
RUN apt-get update && apt-get install -y gcc
RUN head -n 1 /requirements/default.txt | xargs pip install
RUN pip install -r /requirements/default.txt
RUN apt-get remove --purge -y gcc && apt-get autoremove -y
