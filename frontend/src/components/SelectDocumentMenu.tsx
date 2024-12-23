import { createHideableComponent } from "@react-aria/collections";
import useAppState from "../appState";
import { Select, SelectItem } from "./Select";

export const SelectDocumentMenu = createHideableComponent(() => {
  const { appState, setAppState } = useAppState();
  return (
    <Select
      id="select-file"
      wide
      selectedKey={appState.activeDocumentId}
      onSelectionChange={(value) =>
        setAppState((prev) => ({
          ...prev,
          activeDocumentId: value as string,
        }))
      }
      items={appState.documents.map(({ id, basename }) => ({
        id,
        name: basename,
      }))}
    >
      {(item) => <SelectItem textValue={item.name}>{item.name}</SelectItem>}
    </Select>
  );
});
