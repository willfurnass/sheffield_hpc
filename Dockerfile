FROM plaindocs/docker-sphinx

RUN apt-get update && apt-get install -y -q git python-pip

RUN pip install sphinx_bootstrap_theme

RUN mkdir /root/.ssh
ADD drone_ssh_config /root/.ssh/config
ADD drone_id_rsa /root/.ssh/id_rsa
ADD drone_id_rsa.pub /root/.ssh/id_rsa.pub
ADD sync_built_docs.sh /root/sync_built_docs.sh


CMD ["/bin/bash"]
