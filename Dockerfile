FROM python:3.9
RUN apt update && apt install -y libgeos-dev
WORKDIR /src/link-budget/
ADD . .
RUN make && make install
RUN pip3 install -r test_requirements.txt
ENTRYPOINT ["bash"]
