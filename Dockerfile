FROM plaindocs/docker-sphinx

RUN apt-get update && apt-get install -y -q git

CMD ["/bin/bash"]
