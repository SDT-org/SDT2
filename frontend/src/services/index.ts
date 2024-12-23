import * as dataServices from "./data";
import * as docServices from "./documents";
import * as fileServices from "./files";
import * as runServices from "./run";

export const services = {
  ...fileServices,
  ...docServices,
  ...dataServices,
  ...runServices,
};
