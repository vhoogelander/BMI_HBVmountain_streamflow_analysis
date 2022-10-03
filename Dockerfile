FROM continuumio/anaconda3
MAINTAINER Vincent Hoogelander <v.hoogelander@student.tudelft.nl>

# Install here your BMI model:
RUN git clone https://github.com/vhoogelander/Thesis.git /opt/HBVmountain
WORKDIR /opt/HBVmountain

RUN conda env create -f environment_HBVmountain.yml

ENTRYPOINT ["conda", "run", "-n", "HBVmountain", "python3", "-m", "Container.BMI_HBVmountain_Python"]