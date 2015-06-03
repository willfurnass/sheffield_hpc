FROM plaindocs/docker-sphinx

RUN apt-get update && apt-get install -y -q git python-pip

RUN pip install sphinx_bootstrap_theme

CMD ["/bin/bash"]
