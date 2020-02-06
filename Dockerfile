FROM alpine:latest

RUN apk add python3 py3-pip build-base python3-dev zlib-dev && python3 -m pip install pybigwig numpy && apk del py3-pip build-base
COPY src/app/ /app
COPY src/scripts/* /bin
