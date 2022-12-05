FROM continuumio/anaconda3

MAINTAINER Vincent Hoogelander <v.hoogelander@student.tudelft.nl>

# Install grpc4bmi
RUN pip install grpc4bmi==0.2.3

# Install BMI model:
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
RUN pip install rasterio
RUN pip install shapely
RUN pip install geopandas
RUN pip install suntime
RUN pip install cma
RUN pip install ruamel.yaml


#install julia packages
WORKDIR /opt/HBVmountain/Container/Refactoring
RUN julia install.jl


WORKDIR /opt/HBVmountain/Container

##install gcc compiler


RUN apt-get update && \    
    apt-get -y install gcc mono-mcs && \
    rm -rf /var/lib/apt/lists/*dod
RUN apt-get update && \
    apt-get -y install gcc mono-mcs && \
    rm -rf /var/lib/apt/lists/*

##install custom system image
RUN python3 -m julia.sysimage sys.so
RUN julia-py --sysimage sys.so
RUN pip install protobuf==3.20.*

WORKDIR /opt/HBVmountain/Container
# Run bmi server
ENV BMI_MODULE=BMI_HBVmountain_Python
ENV BMI_CLASS=BMI_HBVmountain
ENV BMI_PORT=55555
#ENTRYPOINT ["run-bmi-server", "--name", "BMI_HBVmountain.BMI_HBVmountain_Python", "--path", /opt/HBVmountain/Container]

ENTRYPOINT ["run-bmi-server", "--path", "/opt/HBVmountain/Container/"]


# Expose the magic grpc4bmi port
EXPOSE 55555

