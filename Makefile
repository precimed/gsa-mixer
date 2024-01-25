all: mixer.sif

mixer.sif: Dockerfile
	docker build -t mixer -f Dockerfile . && scripts/from_docker_image.sh mixer

