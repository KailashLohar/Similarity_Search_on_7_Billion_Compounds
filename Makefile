IMAGE_NAME := kailash-spark-image
CONTAINER_NAME := ${USER}-spark-jupyter
HOST_PORT := $(shell shuf -i 8000-17999 -n 1)

PG_IMAGE_NAME := kailash-pg-image
PG_CONTAINER_NAME := ${USER}-pg-jupyter
PG_PORT := 5432

SL_IMAGE_NAME := kailash-sl-image
SL_CONTAINER_NAME := ${USER}-sl-jupyter
SL_PORT := 8501

SERVICE_NAME := ${USER}-spark-jupyter

build-image:
	@cd devops && \
		IMAGE_NAME=$(IMAGE_NAME) CONTAINER_NAME=$(CONTAINER_NAME) HOST_PORT=$(HOST_PORT) \
		PG_IMAGE_NAME=$(PG_IMAGE_NAME) PG_CONTAINER_NAME=$(PG_CONTAINER_NAME) PG_PORT=$(PG_PORT) \
		SL_IMAGE_NAME=$(SL_IMAGE_NAME) SL_CONTAINER_NAME=$(SL_CONTAINER_NAME) SL_PORT=$(SL_PORT) \
		docker-compose build 

start-container:
	@docker ps --format '{{.Names}}' | grep -q "^$(CONTAINER_NAME)$$" && echo "Already running container: \033[1;32m$(CONTAINER_NAME)\033[0m" || true
	@docker ps --format '{{.Names}}' | grep -q "^$(PG_CONTAINER_NAME)$$" && echo "Already running container: \033[1;32m$(PG_CONTAINER_NAME)\033[0m" || true
	@docker ps --format '{{.Names}}' | grep -q "^$(SL_CONTAINER_NAME)$$" && echo "Already running container: \033[1;32m$(SL_CONTAINER_NAME)\033[0m" || true
	@if ! docker ps --format '{{.Names}}' | grep -q -e "^$(CONTAINER_NAME)$$" -e "^$(PG_CONTAINER_NAME)$$" -e "^$(SL_CONTAINER_NAME)$$"; then \
		cd devops && \
			IMAGE_NAME=$(IMAGE_NAME) CONTAINER_NAME=$(CONTAINER_NAME) HOST_PORT=$(HOST_PORT) \
			PG_IMAGE_NAME=$(PG_IMAGE_NAME) PG_CONTAINER_NAME=$(PG_CONTAINER_NAME) PG_PORT=$(PG_PORT) \
			SL_IMAGE_NAME=$(SL_IMAGE_NAME) SL_CONTAINER_NAME=$(SL_CONTAINER_NAME) SL_PORT=$(SL_PORT) \
			docker compose -p $(SERVICE_NAME) up -d > /dev/null 2>&1 && \
		echo "Successfully started container: \033[1;32m$(CONTAINER_NAME)\033[0m"; \
		sleep 1; \
		URL="http://127.0.0.1:$(HOST_PORT)"; \
		echo "JupyterLab is running at: \033[1;34m$$URL\033[0m"; \
		echo "Successfully started container: \033[1;32m$(PG_CONTAINER_NAME)\033[0m"; \
		echo "Successfully started container: \033[1;32m$(SL_CONTAINER_NAME)\033[0m"; \
		sleep 1; \
		URL_st="http://127.0.0.1:$(SL_PORT)"; \
		echo "Streamlit is running at: \033[1;34m$$URL_st\033[0m"; \
	fi

enter-spark-container:
	@echo "You are inside the Container: \033[1;33m$(CONTAINER_NAME)\033[0m"
	@docker exec -u root -it $(CONTAINER_NAME) bash || true

enter-pg-container:
	@echo "You are inside the Container: \033[1;33m$(PG_CONTAINER_NAME)\033[0m"
	@docker exec -u root -it $(PG_CONTAINER_NAME) bash || true

enter-sl-container:
	@echo "You are inside the Container: \033[1;33m$(SL_CONTAINER_NAME)\033[0m"
	@docker exec -u root -it $(SL_CONTAINER_NAME) bash || true

stop-container:
	@if docker ps -a --format '{{.Names}}' | grep -q -e "^$(CONTAINER_NAME)$$" -e "^$(PG_CONTAINER_NAME)$$" -e "^$(SL_CONTAINER_NAME)$$"; then \
		echo "Stopped and removed container: \033[1;31m$(CONTAINER_NAME)\033[0m"; \
		docker rm -f $(CONTAINER_NAME) > /dev/null 2>&1; \
		echo "Stopped and removed container: \033[1;31m$(PG_CONTAINER_NAME)\033[0m"; \
		docker rm -f $(PG_CONTAINER_NAME) > /dev/null 2>&1; \
		echo "Stopped and removed container: \033[1;31m$(SL_CONTAINER_NAME)\033[0m"; \
		docker rm -f $(SL_CONTAINER_NAME) > /dev/null 2>&1; \
	else echo "\033[1;31mThere are no running containers to stop.\033[0m"; fi
