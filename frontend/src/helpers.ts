export const assertDefined = <T>(value: T) => {
  if (value === undefined) {
    throw new Error(`Expected ${value} to be defined`);
  } else {
    return value;
  }
};
