import * as dataServices from "./data";
import * as fileDialogServices from "./fileDialog";
export const services = {
  ...fileDialogServices,
  ...dataServices,
};
