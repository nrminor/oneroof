//
// This file holds a couple Groovy functions for counting amplicons and reference contigs
//

import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

class Utils {

    public static Integer countFastaHeaders(String filePath) {
        def headerCount = 0
        def file = new File(filePath)

        if (!file.exists()) {
            throw new FileNotFoundException("The file $filePath does not exist.")
        }

        file.eachLine { line ->
            if (line.trim().startsWith('>')) {
                headerCount++
            }
        }

        return headerCount
    }

    public static Integer countAmplicons(String filePath) {
        // Read all lines from the file
        def lines = Files.readAllLines(Paths.get(filePath))

        // Extract the fourth column, replace hyphens, and remove _LEFT and _RIGHT
        def processedItems = lines.collect { line ->
            def columns = line.split('\t')
            if (columns.size() >= 4) {
                def item = columns[3]
                item = item.replace('-', '_')
                item = item.replace('_LEFT', '').replace('_RIGHT', '')
                return item
            }
            return null
        }.findAll { it != null }

        // Convert to a set to find unique items
        def uniqueItems = processedItems.toSet()

        // Return the number of unique items
        return uniqueItems.size()
        }

    }
