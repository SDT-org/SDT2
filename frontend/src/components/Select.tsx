import type {
  ListBoxItemProps,
  SelectProps as RACSelectProps,
  ValidationResult,
} from "react-aria-components";
import {
  Button,
  FieldError,
  Label,
  ListBox,
  ListBoxItem,
  Popover,
  Select as RACSelect,
  SelectValue,
  Text,
} from "react-aria-components";

import type React from "react";

interface SelectProps<T extends object>
  extends Omit<RACSelectProps<T>, "children"> {
  label?: string;
  description?: string;
  errorMessage?: string | ((validation: ValidationResult) => string);
  items?: Iterable<T>;
  children: React.ReactNode | ((item: T) => React.ReactNode);
}

export function Select<T extends object>({
  label,
  description,
  errorMessage,
  children,
  items,
  ...props
}: SelectProps<T>) {
  return (
    <RACSelect {...props}>
      <Label>{label}</Label>
      <Button>
        <SelectValue />
        <svg
          className="caret"
          aria-hidden="true"
          xmlns="http://www.w3.org/2000/svg"
          width="16"
          height="16"
          fill="none"
          viewBox="0 0 24 24"
        >
          <path
            stroke="currentColor"
            stroke-linecap="round"
            stroke-linejoin="round"
            stroke-width="2"
            d="m8 10 4 4 4-4"
          />
        </svg>
      </Button>
      {description && <Text slot="description">{description}</Text>}
      <FieldError>{errorMessage}</FieldError>
      <Popover>
        {items ? <ListBox items={items}>{children}</ListBox> : null}
      </Popover>
    </RACSelect>
  );
}

export function SelectItem(props: ListBoxItemProps) {
  return (
    <ListBoxItem
      {...props}
      className={({ isFocused, isSelected }) =>
        `my-item ${isFocused ? "focused" : ""} ${isSelected ? "selected" : ""}`
      }
    />
  );
}
