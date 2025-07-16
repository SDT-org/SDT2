#!/usr/bin/env bun
import { $ } from "bun";
import {
  existsSync,
  mkdirSync,
  rmSync,
  copyFileSync,
  chmodSync,
} from "node:fs";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import { platform } from "node:os";
import { spawn } from "node:child_process";

const __dirname = dirname(fileURLToPath(import.meta.url));
const isWindows = platform() === "win32";
const isMac = platform() === "darwin";
const isLinux = platform() === "linux";

// Build configuration
const config = {
  lzaniRepo: "https://github.com/omaralvarez/lz-ani.git", // Update with actual repo
  lzaniVersion: "main",
  msys2Packages: [
    "mingw-w64-x86_64-gcc",
    "mingw-w64-x86_64-make",
    "mingw-w64-x86_64-cmake",
    "mingw-w64-x86_64-zlib",
    "git",
    "make",
  ],
  requiredDlls: ["libwinpthread-1.dll", "zlib1.dll"],
  paths: {
    build: join(__dirname, "build"),
    backendBin: join(__dirname, "backend", "bin"),
    temp: join(__dirname, "build", "temp"),
    lzaniSrc: join(__dirname, "build", "temp", "lz-ani"),
    venv: join(__dirname, "venv"),
  },
};

// Utility functions
const ensureDir = (path) => {
  if (!existsSync(path)) {
    mkdirSync(path, { recursive: true });
  }
};

const cleanBuild = () => {
  console.log("🧹 Cleaning build artifacts...");
  if (existsSync(config.paths.build)) {
    rmSync(config.paths.build, { recursive: true, force: true });
  }
  console.log("✅ Clean completed");
};

const checkMsys2 = () => {
  if (!isWindows) return null;

  const msys2Paths = [
    "C:\\msys64",
    "C:\\msys2",
    `${process.env.PROGRAMFILES}\\msys64`,
    `${process.env["PROGRAMFILES(X86)"]}\\msys64`,
  ];

  for (const path of msys2Paths) {
    if (existsSync(path)) {
      console.log(`✅ Found MSYS2 at: ${path}`);
      return path;
    }
  }

  throw new Error(
    "❌ MSYS2 not found. Please install from https://www.msys2.org/",
  );
};

const runMsys2Command = async (command, msys2Path) => {
  if (!isWindows) {
    return $`sh -c ${command}`;
  }

  const msys2Shell = join(msys2Path, "msys2_shell.cmd");
  return new Promise((resolve, reject) => {
    const proc = spawn(
      msys2Shell,
      ["-defterm", "-no-start", "-mingw64", "-c", command],
      {
        stdio: "inherit",
      },
    );
    proc.on("close", (code) => {
      if (code === 0) resolve();
      else reject(new Error(`Command failed with code ${code}`));
    });
  });
};

const installMsys2Packages = async (msys2Path) => {
  if (!isWindows) return;

  console.log("📦 Installing MSYS2 packages...");
  const packages = config.msys2Packages.join(" ");
  await runMsys2Command(
    `pacman -S --noconfirm --needed ${packages}`,
    msys2Path,
  );
};

const cloneLzaniSource = async () => {
  ensureDir(config.paths.temp);

  if (existsSync(config.paths.lzaniSrc)) {
    console.log("🔄 Updating lz-ani source...");
    await $`cd ${config.paths.lzaniSrc} && git pull`;
  } else {
    console.log("📥 Cloning lz-ani source...");
    await $`cd ${config.paths.temp} && git clone ${config.lzaniRepo} lz-ani`;
  }

  await $`cd ${config.paths.lzaniSrc} && git checkout ${config.lzaniVersion}`;
};

const createLzaniMakefile = async () => {
  const makefilePath = join(config.paths.lzaniSrc, "Makefile");
  if (existsSync(makefilePath)) return;

  console.log("📝 Creating Makefile for lz-ani...");
  const makefile = `
# Makefile for lz-ani
CC = gcc
CXX = g++
CFLAGS = -Wall -O3 -std=c11 -fopenmp
CXXFLAGS = -Wall -O3 -std=c++17 -fopenmp
LDFLAGS = -lz -lm -fopenmp

# Windows specific flags
ifeq ($(OS),Windows_NT)
    LDFLAGS += -static-libgcc -static-libstdc++
    EXE_EXT = .exe
else
    EXE_EXT =
endif

# Source files
SOURCES = $(wildcard src/*.c src/*.cpp)
OBJECTS = $(SOURCES:.c=.o)
OBJECTS := $(OBJECTS:.cpp=.o)

# Target
TARGET = lz-ani$(EXE_EXT)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
`;

  await Bun.write(makefilePath, makefile);
};

const buildLzani = async (msys2Path) => {
  console.log("🔨 Building lz-ani...");
  await createLzaniMakefile();

  if (isWindows) {
    await runMsys2Command(
      `cd '${config.paths.lzaniSrc}' && make clean && make`,
      msys2Path,
    );
  } else {
    await $`cd ${config.paths.lzaniSrc} && make clean && make -j4`;
  }
};

const copyLzaniBinary = async () => {
  const exeName = isWindows ? "lz-ani.exe" : "lz-ani";
  let srcPath = join(config.paths.lzaniSrc, exeName);

  // Try common build output directories
  if (!existsSync(srcPath)) {
    for (const subdir of ["", "build", "bin", "Release", "Debug"]) {
      const altPath = join(config.paths.lzaniSrc, subdir, exeName);
      if (existsSync(altPath)) {
        srcPath = altPath;
        break;
      }
    }
  }

  if (!existsSync(srcPath)) {
    throw new Error(`Could not find built lz-ani executable at ${srcPath}`);
  }

  const dstPath = join(config.paths.backendBin, exeName);
  console.log(`📋 Copying ${srcPath} to ${dstPath}`);
  ensureDir(config.paths.backendBin);
  copyFileSync(srcPath, dstPath);

  if (!isWindows) {
    chmodSync(dstPath, 0o755);
  }
};

const copyWindowsDlls = async (msys2Path) => {
  if (!isWindows) return;

  console.log("📋 Copying required DLLs...");
  const mingwBin = join(msys2Path, "mingw64", "bin");

  for (const dll of config.requiredDlls) {
    const srcDll = join(mingwBin, dll);
    if (existsSync(srcDll)) {
      const dstDll = join(config.paths.backendBin, dll);
      console.log(`  • ${dll}`);
      copyFileSync(srcDll, dstDll);
    } else {
      console.warn(`⚠️  ${dll} not found in ${mingwBin}`);
    }
  }
};

const buildMainApp = async () => {
  console.log("\n🏗️  Building main SDT2 application...");
  const pythonCmd = isWindows ? "python" : "python3";
  await $`${pythonCmd} build.py`;
};

const createNsisInstaller = async () => {
  if (!isWindows) {
    console.log("ℹ️  NSIS installer is only for Windows");
    return;
  }

  const nsisPaths = [
    "C:\\Program Files\\NSIS\\makensis.exe",
    "C:\\Program Files (x86)\\NSIS\\makensis.exe",
  ];

  let makensis = null;
  for (const path of nsisPaths) {
    if (existsSync(path)) {
      makensis = path;
      break;
    }
  }

  // Try in PATH
  try {
    await $`where makensis`;
    makensis = "makensis";
  } catch {}

  if (!makensis) {
    console.warn("⚠️  NSIS not found. Skipping installer creation.");
    console.log("   Install NSIS from https://nsis.sourceforge.io/");
    return;
  }

  console.log("\n📦 Creating NSIS installer...");
  await $`${makensis} installer.nsis`;
  console.log("✅ Installer created: SDT2_Installer.exe");
};

const checkLzaniBinaries = () => {
  const exeName = isWindows ? "lz-ani.exe" : "lz-ani";
  const lzaniPath = join(config.paths.backendBin, exeName);

  if (!existsSync(lzaniPath)) {
    console.warn(`⚠️  ${exeName} not found in backend/bin/`);
    return false;
  }

  if (isWindows) {
    for (const dll of config.requiredDlls) {
      const dllPath = join(config.paths.backendBin, dll);
      if (!existsSync(dllPath)) {
        console.warn(`⚠️  ${dll} not found in backend/bin/`);
      }
    }
  }

  console.log("✅ lz-ani binaries check passed");
  return true;
};

// Build commands
const commands = {
  async all() {
    console.log("🚀 Starting full build...\n");
    await commands.lzani();
    await commands.main();
    await commands.installer();
  },

  async lzani() {
    console.log("🧬 Building lz-ani...\n");
    let msys2Path = null;

    if (isWindows) {
      msys2Path = checkMsys2();
      await installMsys2Packages(msys2Path);
    }

    await cloneLzaniSource();
    await buildLzani(msys2Path);
    await copyLzaniBinary();

    if (isWindows) {
      await copyWindowsDlls(msys2Path);
    }

    console.log("✅ lz-ani build completed");
  },

  async main() {
    if (!checkLzaniBinaries()) {
      throw new Error(
        "Please build or add lz-ani binaries first (run: bun build.js lzani)",
      );
    }
    await buildMainApp();
    console.log("✅ SDT2 build completed");
  },

  async installer() {
    await createNsisInstaller();
  },

  async clean() {
    cleanBuild();
  },

  async check() {
    checkLzaniBinaries();
  },

  async help() {
    console.log(`
🛠️  SDT2 Build System (Bun)

Usage: bun build.js [command]

Commands:
  all        Build everything (lz-ani + SDT2 + installer)
  lzani      Build lz-ani from source
  main       Build SDT2 application only
  installer  Create NSIS installer (Windows only)
  clean      Clean build artifacts
  check      Check if lz-ani binaries are present
  help       Show this help message

Examples:
  bun build.js all          # Full build
  bun build.js main         # Build SDT2 only (skip lz-ani)
  bun build.js clean        # Clean build files

Platform: ${platform()}
    `);
  },
};

// Main
const main = async () => {
  const command = process.argv[2] || "help";

  if (!commands[command]) {
    console.error(`❌ Unknown command: ${command}`);
    await commands.help();
    process.exit(1);
  }

  try {
    await commands[command]();
  } catch (error) {
    console.error(`\n❌ Build failed: ${error.message}`);
    process.exit(1);
  }
};

main();
