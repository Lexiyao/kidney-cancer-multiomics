FROM bioconductor/bioconductor_docker:RELEASE_3_23

# --- System Python for MOFA2 (reticulate target; basilisk stays external) ---
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Pin mofapy2 + MOFA2 SystemRequirements into the system Python.
# --ignore-installed: the Bioconductor base ships Debian-managed numpy/scipy
# (no RECORD file) that pip cannot uninstall; install the pinned versions over
# them into /usr/local (which precedes dist-packages on sys.path).
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir --break-system-packages --ignore-installed -r /tmp/requirements.txt

# --- MOFA2 basilisk/reticulate resolution (design spec section 5) ---
# Point reticulate at the system Python and force basilisk to an external
# system dir so it never downloads conda inside the container.
ENV RETICULATE_PYTHON=/usr/bin/python3
ENV BASILISK_EXTERNAL_DIR=/opt/basilisk
ENV BASILISK_USE_SYSTEM_DIR=1
RUN mkdir -p /opt/basilisk

WORKDIR /project

# Restore the pinned R library from the lockfile before copying sources.
COPY renv.lock renv.lock
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json
COPY .Rprofile .Rprofile
COPY DESCRIPTION DESCRIPTION
RUN R -e "options(renv.config.pak.enabled = FALSE); renv::restore(prompt = FALSE)"

COPY . /project

CMD ["R", "--no-save"]
