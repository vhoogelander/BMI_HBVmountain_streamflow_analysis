FROM continuumio/anaconda3
MAINTAINER Vincent Hoogelander <v.hoogelander@student.tudelft.nl>

# Install here your BMI model:
RUN git clone https://github.com/vhoogelander/Thesis.git /HBVmountain
WORKDIR /HBVmountain

RUN conda env create -f environment_HBVmountain.yml
#RUN conda create --name BMI_HBVmountain_Python env -f environment_HBVmountain.yml
#RUN conda create --name BMI_HBVmountain_Python -f environment_HBVmountain.yml
#WORKDIR ~/opt/HBVmountain

ENTRYPOINT ["conda", "run", "-n", "base", "python3", "-m", "Container.BMI_HBVmountain_Python"]