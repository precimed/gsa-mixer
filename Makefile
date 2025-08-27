all: gsa-mixer.sif

gsa-mixer.sif: Dockerfile
	docker build -t gsa-mixer -f Dockerfile . && scripts/from_docker_image.sh gsa-mixer

