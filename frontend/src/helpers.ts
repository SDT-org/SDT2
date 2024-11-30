export const assertDefined = <T>(value: T) => {
  if (value === undefined) {
    throw new Error(`Expected ${value} to be defined`);
  } else {
    return value;
  }
};

// https://stackoverflow.com/a/18650828
export const formatBytes = (bytes: number, decimals = 2) => {
  if (!+bytes) return "0 Bytes";

  const k = 1024;
  const dm = decimals < 0 ? 0 : decimals;
  const sizes = [
    "Bytes",
    "KiB",
    "MiB",
    "GiB",
    "TiB",
    "PiB",
    "EiB",
    "ZiB",
    "YiB",
  ];

  const i = Math.floor(Math.log(bytes) / Math.log(k));

  return `${Number.parseFloat((bytes / Math.pow(k, i)).toFixed(dm))} ${sizes[i]}`;
};

export const formatTitle = (key: string) => key.replaceAll("_", " ");
