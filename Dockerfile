FROM plaindocs/docker-sphinx

RUN apt-get update && apt-get install -y -q git python-pip

RUN pip install sphinx_bootstrap_theme

RUN mkdir ~/.ssh
ADD drone_id_rsa ~/.ssh/id_rsa
ADD sync_built_docs.sh

CMD ["/bin/bash"]
