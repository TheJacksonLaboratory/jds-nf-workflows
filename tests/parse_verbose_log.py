import codecs
import re

def read_unicode_file(filename):
    process_name = None
    results = {}
    with codecs.open(filename, 'r', encoding='utf-8') as f:
        current_test = None
        for line in f:
            processed_line = line.strip()
            # Remove non-ASCII and ANSI escape codes
            processed_line = re.sub(r'^[^\x20-\x7E]+', '', processed_line)
            processed_line = re.sub(r'\x1b\[[0-9;]*[A-Za-z]', '', processed_line)
            # Capture 'Test Process <name>'
            match_process = re.search(r"Test Process (.+)$", processed_line)
            if match_process:
                process_name = match_process.group(1)
                results[process_name] = []
                current_test = None
            # Capture 'Test Process <name>'
            match_process = re.search(r"Test Workflow (.+)$", processed_line)
            if match_process:
                process_name = match_process.group(1)
                results[process_name] = []
                current_test = None
            # Capture 'Test [hash] '<test name>''
            match_test = re.search(r"Test \[[^\]]+\] '([^']+)'", processed_line)
            if match_test and process_name:
                test_name = match_test.group(1)
                results[process_name].append({'test': test_name, 'warns': [], 'status': None})
                current_test = results[process_name][-1]
            # Capture WARN for undefined parameter
            match_warn = re.search(r"WARN: Access to undefined parameter `([^`]+)`", processed_line)
            if match_warn and current_test is not None:
                param_name = match_warn.group(1)
                current_test['warns'].append(param_name)
            # Capture PASSED or FAILED status
            match_status = re.search(r"(PASSED|FAILED)", processed_line)
            if match_status and current_test is not None and current_test['status'] is None:
                current_test['status'] = match_status.group(1)

    # Print results as a table

    table = []
    for proc, tests in results.items():
        for t in tests:
            warns = ", ".join(t['warns']) if t['warns'] else ""
            table.append([proc, t['test'], t['status'], warns])

    headers = ["Process|Workflow", "Test", "Status", "Undefined Parameters (WARN)"]
    print(tabulate(table, headers=headers, tablefmt="github"))

if __name__ == "__main__":
    import sys
    from tabulate import tabulate
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
    else:
        read_unicode_file(sys.argv[1])