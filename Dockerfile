FROM continuumio/anaconda3
MAINTAINER Vincent Hoogelander <v.hoogelander@student.tudelft.nl>



# Instal here your BMI model:
RUN git clone https://github.com/vhoogelander/Thesis.git /opt/HBVmountain
WORKDIR /opt/HBVmountain

#install julia
ARG julia_version=1.7.3
RUN curl https://raw.githubusercontent.com/JuliaCI/install-julia/master/install-julia.sh | sed -E "s/\bsudo +//g" | bash -s $julia_version

RUN pip install julia
RUN python -c 'from julia import install; install()'

RUN pip install bmi-python
RUN pip install bmipy
RUN pip install netCDF4

#install julia packages
WORKDIR /opt/HBVmountain/Container/Refactoring
RUN julia install.jl


WORKDIR /opt/HBVmountain/Container

#install gcc compiler
RUN apt-get update && \    
    apt-get -y install gcc mono-mcs && \
    rm -rf /var/lib/apt/lists/*dod
#install custom system image
RUN python3 -m julia.sysimage sys.so

