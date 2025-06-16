# Contributing to XspecT

## Introduction
Thank you for your interest in contributing to XspecT! This page provides guidelines for contributing to the project, including how to set up your own development environment, the XspecT architecture, CI/CD, and the process for submitting contributions.

When contributing to XspecT, please follow the following steps to ensure a smooth process:

- **Read the documentation**: Familiarize yourself with the project by reading the [documentation](https://bionf.github.io/XspecT2/), including the [Understanding XspecT](understanding.md) page and the [architecture overview](#architecture-overview).
- **Follow the coding standards**: Adhere to the project's coding standards and best practices. This includes using consistent naming conventions, writing clear and concise code, and documentation. Furthermore, please make sure your changes are algined with the project's [architecture](#architecture-overview).
- **Write tests**: Ensure that your changes are covered by tests. We use [pytest](https://docs.pytest.org/en/stable/) for testing. If you add new features or fix bugs, please include tests to verify your changes.
- **Document your changes**: Update the documentation to reflect any new features or changes you make. This includes updating the README, Google-style docstrings, and the [Mkdocs](https://www.mkdocs.org)-based documentation.
- **Use clear commit messages**: When committing your changes, use clear and descriptive commit messages that explain the purpose of the changes.
- **Follow the pull request process**: When you're ready to submit your changes, follow the [pull request process](#pull-request-process) outlined below.

## Development Installation
To set up XspecT for development, first make sure you have [Python](https://www.python.org/downloads/) and [Node.js](https://nodejs.org/en/download/) installed. Please note that XspecT is currently not supported in Windows or Alpine Linux environments, unless you build [COBS](https://github.com/aromberg/cobs) yourself.

Get started by cloning the repository:
```bash
git clone https://github.com/BIONF/XspecT2.git
```

You then need to build the web application using Vite. Navigate to the `xspect-web` directory and run the build command, which will also watch for changes:
```bash
cd XspecT2/src/xspect/xspect-web
```
```bash
npx vite build --watch
```

Finally, in a separate terminal, navigate to the root of the cloned repository and install the Python package in editable mode:
```bash
pip install -e .
```

By combining the two processes, you can develop both the frontend and backend simultaneously.

## Architecture Overview
XspecT consists of a Python component (`src/xspect`) and a web application built with [Vite](https://vitejs.dev/) (`src/xspect/xspect-web`). The Python component provides the core functionality, including the command-line interface (CLI) and the backend API, while the web application provides a user-friendly interface for interacting with XspecT. Furthermore, tests for the Python component reside in the `tests/` directory, while documentation is provided in the `docs/` directory.

### Python Component

The Python component of XspecT is structured as follows:

- `main.py`: The entry point for the command-line interface (CLI) and the backend API.
- `web.py`: The [FastAPI](https://fastapi.tiangolo.com/) application that serves the web interface and handles API requests.

The core functionality of XspecT is implemented using the following modules:

- `classify.py`: Contains methods to classify sequences based on previously trained XspecT models.
- `filter_sequences.py`: Contains methods to filter sequences based on classification results.
- `model_management.py`: Contains methods to manage XspecT models.
- `train.py`: Contains methods to train XspecT models based on user-provided data or data from the NCBI/PubMLST API.
- `download_models.py`: Contains methods to download pre-trained XspecT models.

In the background, these modules utilize model classes and a result class, which are defined in the `/models/` folder.

- `/models/probabilistic_filter_model.py`: Base class for probabilistic filter models, which uses COBS indices for classification and stores the model's metadata. Results from the classification are stored in a `ModelResult` class.
- `/models/probabilistic_filter_svm_model.py`: This class extends the base model class and implements a probabilistic filter model, in which classification scores are passed to a support vector machine (SVM) for a final prediction. This model is typically used for species-level classification.
- `/models/probabilistic_filter_mlst_model.py`: This class extends the base model class and implements multilocus strain typing (MLST) by using multiple COBS indices.
- `/models/probabilistic_single_filter_model.py`: This class extends the base model class and implements a model that uses a single Bloom filter for classification. It is typically used for genus-level classification.
- `/models/result.py`: Contains the `ModelResult` class, which stores the results of a classification operation, including classification metadata, hits, and a prediction, if applicable.

Supplementary modules are documented in their respective files.

### Web Application
The web application (`src/xspect/xspect-web`) is built using Vite, [Axios](https://axios-http.com/), [Tailwind CSS](https://tailwindcss.com/), and [shadcn/ui](https://ui.shadcn.com/). It provides a user-friendly interface for interacting with XspecT and includes the following main components:

- `src/api.ts`: Contains the API client for making requests to the backend FastAPI application.
- `src/App.tsx`: The main application component that renders the user interface. It uses React Router for navigation and includes the main layout as well as routing logic.
- `src/assets/`: Contains static assets such as images and icons.
- `src/components/`: Contains reusable components for the user interface, such as buttons, forms, and modals.
- `src/components/ui/`: Contains UI components from shadcn/ui, which are used to build the user interface.
- `src/types.ts`: Contains TypeScript type definitions for the application, including types for API responses.
- `vite.config.ts`: The Vite configuration file that defines how the web application is built and served. Also includes a configuration for the API proxy to the FastAPI backend.

## Continuous Integration and Deployment
We use GitHub Actions to run checks on commits and pull requests. These checks include:

- **Code style and formatting**: Ensures that changes align with the project's code style. We use [Black](https://black.readthedocs.io/en/stable/) for Python code formatting.
- **Linting**: [Pylint](https://pylint.pycqa.org/en/latest/) is used for Python code linting. It checks for coding standards, potential errors, and code smells.
- **Tests**: Ensures that all tests pass. We use [pytest](https://docs.pytest.org/en/stable/) for testing.

Additionally, Github Actions are also used for deployment:

- **Documentation**: The Mkdocs-based documentation is built and deployed to GitHub Pages on changes to the `main` branch. You can view the documentation at [https://bionf.github.io/XspecT2/](https://bionf.github.io/XspecT2/).
- **Python package**: The Python package is built and uploaded to PyPI when a new release is created. This allows users to easily install the latest version of XspecT using `pip install xspect`. Pre-releases are uploaded to TestPyPI and can be installed using `pip install --index-url https://test.pypi.org/simple/ xspect`.

## Pull Request Process
Once you have made your changes and tested them, you can submit a pull request. Please follow these steps:

1. Ensure your code is up to date with the `dev` branch
2. Create a pull request with a clear description of your changes to the `dev` branch
3. Address any feedback from reviewers
4. Once approved, your changes will be merged