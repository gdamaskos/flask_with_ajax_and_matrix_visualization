FROM python:3.8

# install dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -yq apt-utils
RUN DEBIAN_FRONTEND=noninteractive apt-get install -yq htop
RUN apt-get clean
RUN apt-get update && apt-get -qq -y install curl
RUN apt-get -y install ncoils
RUN apt-get -y install muscle
RUN apt-get -y install git
RUN apt -y update
RUN apt -y upgrade
RUN apt -y install python3-pip
RUN apt -y install build-essential libssl-dev libffi-dev python3-dev

# install ccd
RUN pip install --no-cache-dir wheel
WORKDIR /usr/src/flask_with_ajax_and_matrix_visualization
COPY . .
RUN pip install --no-cache-dir -r requirements.txt


# create settings file
RUN touch ./config/ccd_settings.cfg
RUN echo "NCOILS = '/usr/bin/ncoils'" >> ./config/ccd_settings.cfg
RUN echo "IUPRED = '/usr/src/flask_with_ajax_and_matrix_visualization/deps/iupred/iupred'" >> ./config/ccd_settings.cfg
RUN echo "PREDATOR = '/usr/src/flask_with_ajax_and_matrix_visualization/deps/predator/predator'" >> ./config/ccd_settings.cfg
RUN echo "MUSCLE = '/usr/bin/muscle'" >> ./config/ccd_settings.cfg
RUN echo "SECRET_KEY = '\xbb\xca\xca\x9f%:\xca\x17\x03\xf0\xbd\xc8\xbat\xe5\xd2\xb7\x02\xfa\x7f\x8a;\xfc\xf1\xa2\xe4\x02\xcfgo\xd5\xb9'" >> ./config/ccd_settings.cfg
RUN echo "EMAIL = 'proteinccd@gmail.com'" >> ./config/ccd_settings.cfg
RUN echo "BLAST= '/usr/bin/blastp'" >> ./config/ccd_settings.cfg

# define the port number the container should expose
EXPOSE 5000

# final command
CMD ["./run_autoreload.sh"]
