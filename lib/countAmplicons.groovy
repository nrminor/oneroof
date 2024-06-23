import java.nio.file.Files
import java.nio.file.Paths

class AmpliconCounter {
    /**
     * Processes a tab-delimited text file, extracts the fourth column, and returns the number of unique items.
     *
     * @param filePath The path to the tab-delimited text file.
     * @return The number of unique items in the fourth column after processing.
     */
    public static int countFromBed(String filePath) {
        List<String> lines = Files.readAllLines(Paths.get(filePath))

        Set<String> uniqueItems = lines.collect { String line ->
            String[] columns = line.split('\t')
            if (columns.length >= 4) {
                String item = columns[3]
                item = item.replace('-', '_')
                item = item.replace('_LEFT', '').replace('_RIGHT', '')
                return item
            }
            return null
        }.findAll { it != null }.toSet()

        return uniqueItems.size()
    }
}
