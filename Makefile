all: gsa-mixer.sif

gsa-mixer.sif: Dockerfile
	docker build -t gsa-mixer -f Dockerfile . && docker/scripts/from_docker_image.sh gsa-mixer

