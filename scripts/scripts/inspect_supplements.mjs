import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const paths = process.argv.slice(2);

for (const path of paths) {
  const input = await FileBlob.load(path);
  const workbook = await SpreadsheetFile.importXlsx(input);
  const sheets = workbook.worksheets.items;

  console.log(`FILE\t${path}`);
  console.log(`SHEETS\t${sheets.map((sheet) => sheet.name).join("\t")}`);

  for (const sheet of sheets) {
    const used = sheet.getUsedRange(true);
    if (!used) {
      console.log(`EMPTY\t${sheet.name}`);
      continue;
    }

    const rowCount = Math.min(used.rowCount, 100);
    const columnCount = Math.min(used.columnCount, 20);
    const preview = sheet.getRangeByIndexes(
      used.rowIndex,
      used.columnIndex,
      rowCount,
      columnCount,
    );
    console.log(
      `RANGE\t${sheet.name}\t${used.address}\t${used.rowCount}\t${used.columnCount}`,
    );
    console.log(JSON.stringify(preview.values));
  }
}
