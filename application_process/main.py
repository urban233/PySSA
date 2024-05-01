import time
import application_process_manager


def main():
    """Main function."""
    process_manager = application_process_manager.ApplicationProcessManager()
    process_manager.start_pymol()
    while True:
        if process_manager.pymol_process.poll() is not None:
            print("PyMOL crashed!")
            process_manager.start_pymol()
        time.sleep(2)


if __name__ == '__main__':
    main()
