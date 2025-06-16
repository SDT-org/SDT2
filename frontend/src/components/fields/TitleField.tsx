import {
  Input,
  Label,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
} from "react-aria-components";
import type { DocState } from "../../appState";

export const TitleField = ({
  textId,
  textValue,
  onTextChange,
  fontValue,
  onFontChange,
}: {
  textId: string;
  textValue: string;
  onTextChange: (value: string) => void;
  fontValue: string;
  onFontChange: (value: DocState["heatmap"]["titleFont"]) => void;
}) => (
  <>
    <div className="field">
      <Label htmlFor={textId}>Text</Label>
      <TextField id={textId} onChange={onTextChange} value={textValue}>
        <Input />
      </TextField>
    </div>

    <Label htmlFor="font">Font</Label>
    <ToggleButtonGroup
      data-compact
      selectionMode="single"
      disallowEmptySelection={true}
      selectedKeys={[fontValue]}
      onSelectionChange={(selection) =>
        onFontChange(
          selection.values().next().value as DocState["heatmap"]["titleFont"],
        )
      }
    >
      <ToggleButton id={"Sans Serif"}>Sans-serif</ToggleButton>
      <ToggleButton id={"Monospace"}>Monospace</ToggleButton>
    </ToggleButtonGroup>
  </>
);
