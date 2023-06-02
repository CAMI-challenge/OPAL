FROM python:3.11.3-slim

ADD *.py /usr/local/bin/
ADD cami_opal /usr/local/bin/cami_opal
ADD cami_opal/utils /usr/local/bin/cami_opal/utils
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt
