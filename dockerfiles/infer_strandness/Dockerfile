FROM libddocker/r_3.6.1_bioc
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz
RUN tar xzvf kallisto_linux-v0.46.1.tar.gz
RUN chmod -R 755 kallisto
RUN mv kallisto/kallisto /usr/bin/
RUN rm kallisto_linux-v0.46.1.tar.gz
