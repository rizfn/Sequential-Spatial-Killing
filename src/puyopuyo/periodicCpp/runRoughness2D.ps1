# Define constants for steps and L
$steps = 2048
$L_values = @(2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 64, 128)

# Define the range of N_species values
$N_species_values = @(4, 5, 6, 7, 10, 15, 32, 64)

# Define the directory containing the executable
$scriptDir = "C:\GitHub\Sequential-Spatial-Killing\src\puyopuyo\cpp"

# Get the total number of logical processors and define the maximum number of concurrent jobs
$totalProcessors = (Get-WmiObject -Class Win32_ComputerSystem).NumberOfLogicalProcessors
$maxJobs = [math]::Max(1, $totalProcessors - 4)  # Leave 4 processors free

# Loop through the L and N_species values
foreach ($L in $L_values) {
    foreach ($N_species in $N_species_values) {
        # Wait until there is a free slot for another job
        while ((Get-Job | Where-Object { $_.State -eq 'Running' }).Count -ge $maxJobs) {
            Start-Sleep -Seconds 1
        }

        # Start a new job for each combination
        Start-Job -ScriptBlock {
            param($scriptDir, $L, $N_species, $steps)

            # Run the executable with the given parameters
            & "$scriptDir\roughness2D.exe" $L $N_species $steps

            # Check if the process succeeded
            if ($LASTEXITCODE -ne 0) {
                Write-Error "Execution failed for L=$L N_species=$N_species"
            }
        } -ArgumentList $scriptDir, $L, $N_species, $steps
    }
}

# Wait for all jobs to complete
Get-Job | Wait-Job

# Clean up the job list
Get-Job | Remove-Job

Write-Host "All simulations completed."