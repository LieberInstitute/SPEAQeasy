FROM libddocker/ubuntu16.04_base:latest

MAINTAINER Jacob Leonard "leonard.jacob09@gmail.com"

######## Update/Install Core ######## 
RUN apt-get update
RUN pip install numpy
RUN pip install -Iv setuptools==36.0.1
RUN pip install -Iv RSeQC==2.6.4
