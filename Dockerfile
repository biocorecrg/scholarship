FROM rocker/r-ver:3.5

# File Author / Maintainer
MAINTAINER Oliver Gonzalez <ogonzalez@crg.eu>

ARG SHINYPORT=7165

# Install external dependencies 
RUN apt-get update -qq && apt-get install -y --no-install-recommends python curl libcurl4-openssl-dev libssl-dev libsqlite3-dev libxml2-dev qpdf git

RUN mkdir /usr/local/shiny
WORKDIR /usr/local/shiny

RUN mkdir Diagrams

COPY Diagrams /usr/local/shiny/Diagrams

COPY deps.R /usr/local/shiny

RUN Rscript /usr/local/shiny/deps.R


VOLUME /usr/local/shiny/Diagrams/shiny_bookmarks

# Clean cache
RUN apt-get clean
RUN set -x; rm -rf /var/lib/apt/lists/*

EXPOSE 7165


CMD cd /usr/local/shiny/; /usr/local/bin/R -e 'shiny::runApp("Diagrams", host="0.0.0.0", port=7165)'







