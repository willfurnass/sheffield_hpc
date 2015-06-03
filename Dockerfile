FROM plaindocs/docker-sphinx

RUN apt-get install git

CMD ["/bin/bash"]
