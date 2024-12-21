import * as dataServices from "./data";
import * as fileDialogServices from "./fileDialog";
import * as runServices from "./run";
export const services = {
  ...fileDialogServices,
  ...dataServices,
  ...runServices,
};
