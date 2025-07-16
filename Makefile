# Makefile for SDT2 with lz-ani support
# Cross-platform build system

# Detect OS
ifeq ($(OS),Windows_NT)
    PLATFORM := windows
    EXE_EXT := .exe
    PYTHON := python
    RM := del /F /Q
    RMDIR := rmdir /S /Q
    MKDIR := mkdir
    SEP := \\
    INSTALLER := makensis installer.nsis
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        PLATFORM := linux
    endif
    ifeq ($(UNAME_S),Darwin)
        PLATFORM := macos
    endif
    EXE_EXT :=
    PYTHON := python3
    RM := rm -f
    RMDIR := rm -rf
    MKDIR := mkdir -p
    SEP := /
    INSTALLER := echo "NSIS installer is Windows only"
endif

# Directories
BUILD_DIR := build
BACKEND_BIN := backend$(SEP)bin
TEMP_DIR := $(BUILD_DIR)$(SEP)temp

# Build targets
.PHONY: all clean install lzani-check main installer help

# Default target
all: lzani-check main installer

# Help target
help:
	@echo "SDT2 Build System"
	@echo "================="
	@echo ""
	@echo "Available targets:"
	@echo "  make all          - Build everything (default)"
	@echo "  make main         - Build SDT2 application only"
	@echo "  make installer    - Create NSIS installer (Windows only)"
	@echo "  make lzani-check  - Verify lz-ani binaries are present"
	@echo "  make clean        - Clean build artifacts"
	@echo "  make help         - Show this help message"
	@echo ""
	@echo "Platform detected: $(PLATFORM)"

# Check for lz-ani binaries
lzani-check:
	@echo "Checking for lz-ani binaries..."
ifeq ($(PLATFORM),windows)
	@if not exist "$(BACKEND_BIN)$(SEP)lz-ani.exe" ( \
		echo ERROR: lz-ani.exe not found in $(BACKEND_BIN) && \
		echo Please ensure lz-ani.exe is present before building && \
		exit 1 \
	)
	@if not exist "$(BACKEND_BIN)$(SEP)libwinpthread-1.dll" ( \
		echo WARNING: libwinpthread-1.dll not found in $(BACKEND_BIN) \
	)
	@if not exist "$(BACKEND_BIN)$(SEP)zlib1.dll" ( \
		echo WARNING: zlib1.dll not found in $(BACKEND_BIN) \
	)
else
	@if [ ! -f "$(BACKEND_BIN)/lz-ani" ]; then \
		echo "ERROR: lz-ani not found in $(BACKEND_BIN)"; \
		echo "Please ensure lz-ani is present before building"; \
		exit 1; \
	fi
endif
	@echo "lz-ani binaries check passed"

# Build main SDT2 application
main: lzani-check
	@echo "Building SDT2 application..."
	$(PYTHON) build.py
	@echo "SDT2 build completed"

# Create installer (Windows only)
installer:
ifeq ($(PLATFORM),windows)
	@echo "Creating NSIS installer..."
	@if exist "$(BUILD_DIR)$(SEP)app.dist" ( \
		$(INSTALLER) \
	) else ( \
		echo ERROR: Build directory not found. Run 'make main' first. \
	)
else
	@echo "NSIS installer is only available on Windows"
endif

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
ifeq ($(PLATFORM),windows)
	@if exist "$(BUILD_DIR)" $(RMDIR) "$(BUILD_DIR)"
	@if exist "__pycache__" $(RMDIR) "__pycache__"
	@if exist "*.spec" $(RM) *.spec
else
	$(RMDIR) $(BUILD_DIR)
	$(RMDIR) __pycache__
	$(RM) *.spec
endif
	@echo "Clean completed"

# Development targets
.PHONY: dev-setup dev-run

# Setup development environment
dev-setup:
	@echo "Setting up development environment..."
	$(PYTHON) -m venv venv
ifeq ($(PLATFORM),windows)
	@echo "Activate virtual environment with: venv\Scripts\activate"
else
	@echo "Activate virtual environment with: source venv/bin/activate"
endif
	@echo "Then run: pip install -r requirements.txt"

# Run in development mode
dev-run:
	@echo "Running SDT2 in development mode..."
	cd backend && $(PYTHON) src$(SEP)app.py
