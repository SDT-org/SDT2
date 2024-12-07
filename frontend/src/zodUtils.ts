import type { z } from "zod";

// All code from zod_utilz
// https://github.com/JacobWeisenburger/zod_utilz?tab=readme-ov-file#yet-another-library

const omit = <Obj extends Record<string, unknown>>(
  obj: Obj,
  keys: (keyof Obj | string)[],
) => {
  if (keys.length === 0) return obj;
  const entries = Object.entries(obj) as [keyof Obj, unknown][];
  return Object.fromEntries(entries.filter(([key]) => !keys.includes(key)));
};

const pick = <Obj extends Record<string, unknown>>(
  obj: Obj,
  keys: (keyof Obj)[],
) => {
  if (keys.length === 0) return {};
  const entries = Object.entries(obj) as [keyof Obj, unknown][];
  return Object.fromEntries(entries.filter(([key]) => keys.includes(key)));
};

type ObjectIterator<Obj, Result> = (
  value: Obj[keyof Obj],
  key: string,
  collection: Obj,
) => Result;

const mapValues =
  // biome-ignore lint/suspicious/noExplicitAny: library code
    <Obj extends Record<string, any>, Result>(
      fn: ObjectIterator<Obj, Result>,
    ) =>
    (obj?: Obj) => {
      if (!obj) return {};
      const map = Object.keys(obj).reduce((map, key) => {
        map.set(key, fn(obj[key], key, obj));
        return map;
      }, new Map());
      return Object.fromEntries(map);
    };

/**
 * SPR stands for SafeParseResult
 *
 * This enables [optional chaining](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Operators/Optional_chaining) or [nullish coalescing](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Operators/Nullish_coalescing) for `z.SafeParseReturnType`.
 *
 * ### Usage:
 * ```
 * import { zu } from 'zod_utilz'
 * const schema = z.object( { foo: z.string() } )
 * const result = zu.SPR( schema.safeParse( { foo: 42 } ) )
 * const fooDataOrErrors = result.data?.foo ?? result.error?.format().foo?._errors
 * ```
 */
export function SPR<Input, Output>(
  result: z.SafeParseReturnType<Input, Output>,
): {
  success: (typeof result)["success"];
  data: z.SafeParseSuccess<Output>["data"] | undefined;
  error: z.SafeParseError<Input>["error"] | undefined;
} {
  return result.success
    ? { ...result, error: undefined }
    : { ...result, data: undefined };
}

/**
partialSafeParse allows you to get the valid fields even if there was an error in another field

@example
import { zu } from 'zod_utilz'
const userSchema = z.object( { name: z.string(), age: z.number() } )
const result = zu.partialSafeParse( userSchema, { name: null, age: 42 } )
// {
//     successType: 'partial',
//     validData: { age: 42 },
//     invalidData: { name: null },
// }
result.error?.flatten().fieldErrors
// {
//     name: [ 'Expected string, received null' ],
// }
*/
export function partialSafeParse<Schema extends z.AnyZodObject>(
  schema: Schema,
  input: unknown,
): ReturnType<typeof SPR> & {
  successType: "full" | "partial" | "none";
  validData: Partial<z.infer<Schema>>;
  invalidData: Partial<z.infer<Schema>>;
} {
  const result = SPR(schema.safeParse(input));
  if (result.success)
    return {
      ...result,
      successType: "full",
      validData: result.data as Partial<z.infer<Schema>>,
      invalidData: {},
    } as const;

  const { fieldErrors, formErrors } = result.error?.flatten() ?? {};
  if (formErrors?.length)
    return {
      ...result,
      successType: "none",
      validData: {},
      invalidData: {},
    };

  const inputObj = input as z.infer<Schema>;
  const keysWithInvalidData = Object.keys(fieldErrors ?? {});
  const validInput = omit(inputObj, keysWithInvalidData);
  const invalidData = pick(inputObj, keysWithInvalidData) as Partial<
    z.infer<Schema>
  >;

  const validData = schema
    .omit(mapValues(() => true)(fieldErrors))
    .parse(validInput) as Partial<z.infer<Schema>>;

  return {
    ...result,
    successType: "partial",
    validData,
    invalidData,
  };
}
