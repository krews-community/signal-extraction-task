FROM alpine:latest

RUN apk add python3 py3-pip build-base python3-dev zlib-dev git libstdc++ && \
    python3 -m pip install pybigwig numpy git+git://github.com/esnme/ultrajson.git joblib && \
    apk del py3-pip build-base git
COPY src/app/ /app
COPY src/scripts/* /bin/
RUN rm -rf /var/cache/apk/* && cd /

ENV PYTHONNOUSERSITE 1
