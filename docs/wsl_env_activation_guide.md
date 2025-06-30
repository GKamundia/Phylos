# Activating Conda Environments and Python Virtual Environments (venv) in WSL

This guide explains how to activate and use both Conda environments and Python virtual environments (venv) within Windows Subsystem for Linux (WSL), as used in the RVF Nextstrain pipeline setup.

---

## 1. Activating a Conda Environment in WSL

### Prerequisites

- Miniconda or Anaconda must be installed in your WSL environment.
- Your desired environment (e.g., `nextstrain` or `nextstrain-py39`) must already be created.

### Steps

1. **Open your WSL terminal.**
2. **Initialize Conda (if not already done):**
   ```sh
   conda init bash
   exec bash
   ```
   (You only need to do this once after installing Conda.)
3. **Activate your environment:**
   ```sh
   conda activate nextstrain
   # or for a specific Python version:
   conda activate nextstrain-py39
   ```
4. **Verify activation:**
   The prompt should change to show the environment name, e.g., `(nextstrain) user@host:~$`

---

## 2. Activating a Python Virtual Environment (venv) in WSL

### Prerequisites

- Python 3 is installed in your WSL environment.
- A virtual environment has been created (e.g., `/home/youruser/augur-env`).

### Steps

1. **Open your WSL terminal.**
2. **Activate the venv:**
   ```sh
   source /home/youruser/augur-env/bin/activate
   # or, if you are in the venv directory:
   source bin/activate
   ```
3. **Verify activation:**
   The prompt should change to show the venv name, e.g., `(augur-env) user@host:~$`

---

## 3. Switching Between Conda and venv

- You can only have one environment active at a time.
- To switch, deactivate the current environment first:
  ```sh
  conda deactivate
  # or
  deactivate
  ```
- Then activate the other environment as shown above.

---

## 4. Example Usage in the RVF Nextstrain Pipeline

- **Conda:**
  ```sh
  conda activate nextstrain
  augur --version
  ```
- **venv:**
  ```sh
  source /home/youruser/augur-env/bin/activate
  augur --version
  ```

---

## 5. Troubleshooting

- If `conda` is not found, ensure Miniconda/Anaconda is installed and your shell is initialized (`conda init bash`).
- If `activate` does not work, use `source` as shown above.
- If you see conflicting environments, always `deactivate` before switching.

---

**For more details, see the official documentation:**

- [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
- [Python venv Documentation](https://docs.python.org/3/library/venv.html)

/root/miniconda3/bin/conda init bash
exec bash

conda activate nextstrain-py39
