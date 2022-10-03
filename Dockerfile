FROM ubuntu:bionic
MAINTAINER Vincent Hoogelander <v.hoogelander@student.tudelft.nl>

# Install grpc4bmi
RUN pip install git+https://github.com/eWaterCycle/grpc4bmi.git#egg=grpc4bmi

# Install here your BMI model:
RUN git clone https://github.com/vhoogelander/Thesis.git /opt/HBVmountain

# Run bmi server
ENTRYPOINT ["run-bmi-server", "--name", "Thesis.Container.HBVmountain", "--path", "/opt/HBVmountain"]

# Expose the magic grpc4bmi port
EXPOSE 55555
