FROM continuumio/anaconda3
MAINTAINER Vincent Hoogelander <v.hoogelander@student.tudelft.nl>



# Instal here your BMI model:
RUN git clone https://github.com/vhoogelander/Thesis.git /opt/HBVmountain
WORKDIR /opt/HBVmountain

RUN conda env create -f hbvmountain_environment.yml

#ENTRYPOINT ["conda", "run", "-n", "HBVmountain", "python3", "-m", "Container.BMI_HBVmountain_Python"]
#ENTRYPOINT ["python3", "-m", "Container.BMI_HBVmountain_Python"]