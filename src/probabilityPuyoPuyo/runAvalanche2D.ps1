# Define constants for steps and L
$steps = 1024 * 8
$L = 64

# Define the range of N_species values
$N_species_values = @(6, 7)

# Define the number of simulations
$N_sims = 6

# Define the directory containing the executable
$scriptDir = "C:\GitHub\Sequential-Spatial-Killing\src\probabilityPuyoPuyo"

# Get the total number of logical processors and define the maximum number of concurrent jobs
$totalProcessors = (Get-WmiObject -Class Win32_ComputerSystem).NumberOfLogicalProcessors
$maxJobs = [math]::Max(1, $totalProcessors - 4)  # Leave 4 processors free

# Loop through the N_species values
foreach ($N_species in $N_species_values) {
    # Loop through the number of simulations
    for ($sim = 1; $sim -le $N_sims; $sim++) {
        # Wait until there is a free slot for another job
        while ((Get-Job | Where-Object { $_.State -eq 'Running' }).Count -ge $maxJobs) {
            Start-Sleep -Seconds 1
        }

        # Start a new job for each simulation
        Start-Job -ScriptBlock {
            param($scriptDir, $L, $N_species, $steps)

            # Run the executable with the given parameters
            & "$scriptDir\avalanche2D.exe" $L $N_species $steps

            # Check if the process succeeded
            if ($LASTEXITCODE -ne 0) {
                Write-Error "Execution failed for N_species=$N_species"
            }
        } -ArgumentList $scriptDir, $L, $N_species, $steps
    }
}

# Wait for all jobs to complete
Get-Job | Wait-Job

# Clean up the job list
Get-Job | Remove-Job

Write-Host "All simulations completed."