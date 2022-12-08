help:
	@echo "> Type \"make build\" to build the docker image"
	@echo "> Type \"make run\" to run the program after the build"

build:
	docker build --tag converter:latest .

run:
	docker run -u $(shell id -u) -v $(shell pwd):/files -ti converter:latest
