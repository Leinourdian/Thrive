using System;
using System.Collections.Generic;
using Godot;
using Newtonsoft.Json;
using Array = Godot.Collections.Array;

/// <summary>
///   Membrane for microbes
/// </summary>
public class Membrane : MeshInstance, IComputedMembraneData
{
    /// <summary>
    ///   This must be big enough that no organelle can be at this position.
    ///   TODO: maybe switching to nullable float would be a good alternative?
    /// </summary>
    public const float INVALID_FOUND_ORGANELLE = -999999.0f;

    [Export]
    public ShaderMaterial? MaterialToEdit;

    private static readonly List<Vector2> PreviewMembraneOrganellePositions = new() { new Vector2(0, 0) };

    /// <summary>
    ///   Stores the generated 2-Dimensional membrane. Needed for contains calculations
    /// </summary>
    private readonly List<Vector2> vertices2D = new();

    // Work buffers used when generating membrane data
    private List<Vector2> previousWorkBuffer = new();
    private List<Vector2> nextWorkBuffer = new();

    private float healthFraction = 1.0f;
    private float wigglyNess = 1.0f;
    private float sizeWigglyNessDampeningFactor = 0.22f;
    private float movementWigglyNess = 1.0f;
    private float sizeMovementWigglyNessDampeningFactor = 0.32f;
    private Color tint = Colors.White;
    private float dissolveEffectValue;

    private MembraneType? type;

#pragma warning disable CA2213
    private Texture? albedoTexture;
    private Texture noiseTexture = null!;
#pragma warning restore CA2213

    private string? currentlyLoadedAlbedoTexture;

    private bool dirty = true;
    private bool radiusIsDirty = true;
    private bool convexShapeIsDirty = true;
    private float cachedRadius;
    private Vector3[] cachedConvexShape = null!;

    /// <summary>
    ///   Amount of segments on one side of the above described
    ///   square. The amount of points on the side of the membrane.
    /// </summary>
    private int membraneResolution = Constants.MEMBRANE_RESOLUTION * 1; //*10;

    /// <summary>
    ///   When true the mesh needs to be regenerated and material properties applied
    /// </summary>
    public bool Dirty
    {
        get => dirty;
        set
        {
            if (value)
            {
                radiusIsDirty = true;
                convexShapeIsDirty = true;
            }

            dirty = value;
        }
    }

    /// <summary>
    ///   Organelle positions of the microbe, needs to be set for the membrane to appear
    /// </summary>
    /// <remarks>
    ///   <para>
    ///     The contents in this list should not be modified, a new list should be assigned.
    ///     TODO: change the type here to be a readonly list
    ///   </para>
    /// </remarks>
    public List<Vector2> OrganellePositions { get; set; } = PreviewMembraneOrganellePositions;

    /// <summary>
    ///   Returns a convex shaped 3-Dimensional array of vertices from the generated <see cref="vertices2D"/>.
    /// </summary>
    /// <remarks>
    ///   <para>
    ///     NOTE: This is not the same as the 3D vertices used for the visuals.
    ///   </para>
    /// </remarks>
    [JsonIgnore]
    public Vector3[] ConvexShape
    {
        get
        {
            if (convexShapeIsDirty)
            {
                if (Dirty)
                    Update();

                float height = 0.1f;

                if (Type.CellWall)
                    height = 0.05f;

                cachedConvexShape = new Vector3[vertices2D.Count];
                for (var i = 0; i < vertices2D.Count; ++i)
                {
                    cachedConvexShape[i] = new Vector3(vertices2D[i].x, height / 2, vertices2D[i].y);
                }

                convexShapeIsDirty = false;
            }

            return cachedConvexShape;
        }
    }

    /// <summary>
    ///   The type of the membrane.
    /// </summary>
    /// <exception cref="InvalidOperationException">When trying to read before this is initialized</exception>
    /// <exception cref="ArgumentNullException">If value is attempted to be set to null</exception>
    public MembraneType Type
    {
        get => type ?? throw new InvalidOperationException("Membrane type has not been set yet");
        set
        {
            if (value == null)
                throw new ArgumentNullException();

            if (type == value)
                return;

            type = value;
            dirty = true;
        }
    }

    /// <summary>
    ///   How healthy the cell is, mixes in a damaged texture. Range 0.0f - 1.0f
    /// </summary>
    public float HealthFraction
    {
        get => healthFraction;
        set
        {
            value = value.Clamp(0.0f, 1.0f);
            if (value == HealthFraction)
                return;

            healthFraction = value;
            ApplyHealth();
        }
    }

    /// <summary>
    ///   How much the membrane wiggles. Used values are 0 and 1
    /// </summary>
    public float WigglyNess
    {
        get => wigglyNess;
        set
        {
            wigglyNess = Mathf.Clamp(value, 0.0f, 1.0f);
            ApplyWiggly();
        }
    }

    public float MovementWigglyNess
    {
        get => movementWigglyNess;
        set
        {
            movementWigglyNess = Mathf.Clamp(value, 0.0f, 1.0f);
            ApplyMovementWiggly();
        }
    }

    public Color Tint
    {
        get => tint;
        set
        {
            // Desaturate it here so it looks nicer (could implement as method that
            // could be called i suppose)

            // According to stack overflow HSV and HSB are the same thing
            value.ToHsv(out var hue, out var saturation, out var brightness);

            value = Color.FromHsv(hue, saturation * 0.75f, brightness,
                Mathf.Clamp(value.a, 0.4f - brightness * 0.3f, 1.0f));

            if (tint == value)
                return;

            tint = value;

            // If we already have created a material we need to re-apply it
            ApplyTint();
        }
    }

    /// <summary>
    ///   Quick radius value for the membrane size
    /// </summary>
    public float EncompassingCircleRadius
    {
        get
        {
            if (radiusIsDirty)
            {
                cachedRadius = CalculateEncompassingCircleRadius();
                radiusIsDirty = false;
            }

            return cachedRadius;
        }
    }

    public float DissolveEffectValue
    {
        get => dissolveEffectValue;
        set
        {
            dissolveEffectValue = value;
            ApplyDissolveEffect();
        }
    }

    public override void _Ready()
    {
        type ??= SimulationParameters.Instance.GetMembrane("single");

        if (MaterialToEdit == null)
            GD.PrintErr("MaterialToEdit on Membrane is not set");

        noiseTexture = GD.Load<Texture>("res://assets/textures/dissolve_noise.tres");

        Dirty = true;
    }

    public override void _Process(float delta)
    {
        if (!Dirty)
            return;

        Update();
    }

    /// <summary>
    ///   Sees if the given point is inside the membrane.
    /// </summary>
    /// <remarks>
    ///   <para>
    ///     This is quite an expensive method as this loops all the vertices
    ///   </para>
    /// </remarks>
    public bool Contains(float x, float y)
    {
        bool crosses = false;

        int n = vertices2D.Count;
        for (int i = 0; i < n - 1; i++)
        {
            if ((vertices2D[i].y <= y && y < vertices2D[i + 1].y) ||
                (vertices2D[i + 1].y <= y && y < vertices2D[i].y))
            {
                if (x < (vertices2D[i + 1].x - vertices2D[i].x) *
                    (y - vertices2D[i].y) /
                    (vertices2D[i + 1].y - vertices2D[i].y) +
                    vertices2D[i].x)
                {
                    crosses = !crosses;
                }
            }
        }

        return crosses;
    }

    /// <summary>
    ///   Finds the point on the membrane nearest to the given point.
    /// </summary>
    /// <remarks>
    ///   <para>
    ///     Used for finding out where to put an external organelle.
    ///   </para>
    ///   <para>
    ///     The returned Vector is in world coordinates (x, 0, z) and
    ///     not in internal membrane coordinates (x, y, 0). This is so
    ///     that gameplay code doesn't have to do the conversion
    ///     everywhere this is used.
    ///   </para>
    /// </remarks>
    public Vector3 GetVectorTowardsNearestPointOfMembrane(float x, float y)
    {
        // Calculate now if dirty to make flagella positioning only have to be done once
        // NOTE: that flagella position should only be read once all organelles that are
        // going to be added / removed on this game update are done.
        if (Dirty)
            Update();

        float organelleAngle = Mathf.Atan2(y, x);

        Vector2 closestSoFar = new Vector2(0, 0);
        float angleToClosest = Mathf.Pi * 2;

        foreach (var vertex in vertices2D)
        {
            if (Mathf.Abs(Mathf.Atan2(vertex.y, vertex.x) - organelleAngle) <
                angleToClosest)
            {
                closestSoFar = new Vector2(vertex.x, vertex.y);
                angleToClosest =
                    Mathf.Abs(Mathf.Atan2(vertex.y, vertex.x) - organelleAngle);
            }
        }

        return new Vector3(closestSoFar.x, 0, closestSoFar.y);
    }

    /// <summary>
    ///   Return the position of the closest organelle to the target point if it is less then a certain threshold away.
    /// </summary>
    public Vector2 FindClosestOrganelles(Vector2 target)
    {
        // The distance we want the membrane to be from the organelles squared.
        float closestSoFar = float.MaxValue;// 4;
        Vector2 closest = new Vector2(INVALID_FOUND_ORGANELLE, INVALID_FOUND_ORGANELLE);

        foreach (var pos in OrganellePositions)
        {
            float lenToObject = (target - pos).LengthSquared();

            if (lenToObject < closestSoFar)//lenToObject < 4 && lenToObject < closestSoFar)
            {
                closestSoFar = lenToObject;
                closest = pos;
            }
        }

        return closest;
    }

    public bool MatchesCacheParameters(ICacheableData cacheData)
    {
        if (cacheData is IComputedMembraneData data)
            return this.MembraneDataFieldsEqual(data);

        return false;
    }

    public long ComputeCacheHash()
    {
        return this.ComputeMembraneDataHash();
    }

    /// <summary>
    ///   Decides where the point needs to move based on the position of the closest organelle.
    /// </summary>
    private static Vector2 GetMovement(Vector2 target, Vector2 closestOrganelle)
    {
        float power = Mathf.Pow(2.7f, -(target - closestOrganelle).Length() / 10) / 50;

        return (closestOrganelle - target) * power;
    }

    private static Vector2 GetMovementForCellWall(Vector2 target, Vector2 closestOrganelle)
    {
        float power = Mathf.Pow(10.0f, -(target - closestOrganelle).Length()) / 50;

        return (closestOrganelle - target) * power;
    }

    /// <summary>
    ///   First generates the 2D vertices and then builds the 3D mesh
    /// </summary>
    private void InitializeMesh()
    {
        // First try to get from cache as it's very expensive to generate the membrane
        var cached = this.FetchDataFromCache(ProceduralDataCache.Instance.ReadMembraneData);

        if (cached != null)
        {
            CopyMeshFromCache(cached);
            return;
        }

        var tim1 = Time.GetTicksUsec();

        // The length in pixels (probably not accurate?) of a side of the square that bounds the membrane.
        // Half the side length of the original square that is compressed to make the membrane.
        float cellDimensions = 10; //int to float if you want to mess around

        foreach (var pos in OrganellePositions)
        {
            if (Mathf.Abs(pos.x) + 1 > cellDimensions)
            {
                cellDimensions = (int)Mathf.Abs(pos.x) + 1;
            }

            if (Mathf.Abs(pos.y) + 1 > cellDimensions)
            {
                cellDimensions = (int)Mathf.Abs(pos.y) + 1;
            }
        }

        previousWorkBuffer.Capacity = vertices2D.Capacity;
        nextWorkBuffer.Capacity = previousWorkBuffer.Capacity;

        for (int i = membraneResolution; i > 0; i--)
        {
            previousWorkBuffer.Add(new Vector2(-cellDimensions,
                cellDimensions - 2 * cellDimensions / membraneResolution * i));
        }

        for (int i = membraneResolution; i > 0; i--)
        {
            previousWorkBuffer.Add(new Vector2(
                cellDimensions - 2 * cellDimensions / membraneResolution * i,
                cellDimensions));
        }

        for (int i = membraneResolution; i > 0; i--)
        {
            previousWorkBuffer.Add(new Vector2(cellDimensions,
                -cellDimensions + 2 * cellDimensions / membraneResolution * i));
        }

        for (int i = membraneResolution; i > 0; i--)
        {
            previousWorkBuffer.Add(new Vector2(
                -cellDimensions + 2 * cellDimensions / membraneResolution * i,
                -cellDimensions));
        }

        for (int i = 0; i < previousWorkBuffer.Count; i++)
        {
            previousWorkBuffer[i] = previousWorkBuffer[i] * 100.0f;
        }

        var tim2 = Time.GetTicksUsec();

        DrawCorrectMembrane(cellDimensions, previousWorkBuffer, nextWorkBuffer);
        (previousWorkBuffer, nextWorkBuffer) = (nextWorkBuffer, previousWorkBuffer);

        //// This needs to actually run a bunch of times as the points moving towards the organelles is iterative.
        //// We use rotating work buffers to save time on skipping useless copies
        //for (int i = 0; i < 40 * cellDimensions; i++)// cellDimensions; i++)
        //{
        //    DrawCorrectMembrane(cellDimensions, previousWorkBuffer, nextWorkBuffer);

        //    (previousWorkBuffer, nextWorkBuffer) = (nextWorkBuffer, previousWorkBuffer);
        //}

        var tim3 = Time.GetTicksUsec();

        // Copy final vertex data from the work buffer
        vertices2D.Clear();

        // The work buffer not being pointed to as the next, is the one we should read the result from
        vertices2D.AddRange(previousWorkBuffer);

        previousWorkBuffer.Clear();
        nextWorkBuffer.Clear();

        var tim4 = Time.GetTicksUsec();

        BuildMesh();

        var tim5 = Time.GetTicksUsec();
        //GD.Print("first " + (tim2 - tim1) + ", second " + (tim3 - tim2) + ", third " + (tim4 - tim3) + ", fourth " + (tim5 - tim4));
        GD.Print("========================================================================================" + (tim3 - tim2));
    }

    private int InitializeCorrectMembrane(int writeIndex, Vector3[] vertices,
        Vector2[] uvs)
    {
        // common variables
        float height = 0.1f;
        float multiplier = 2.0f * Mathf.Pi;
        var center = new Vector2(0.5f, 0.5f);

        // cell walls need obvious inner/outer membranes (we can worry
        // about chitin later)
        if (Type.CellWall)
        {
            height = 0.05f;
            multiplier = Mathf.Pi * 2.0f;
        }

        vertices[writeIndex] = new Vector3(0, height / 2, 0);
        uvs[writeIndex] = center;
        ++writeIndex;

        for (int i = 0, end = vertices2D.Count; i < end + 1; i++)
        {
            // Finds the UV coordinates be projecting onto a plane and
            // stretching to fit a circle.

            float currentRadians = multiplier * i / end;

            vertices[writeIndex] = new Vector3(vertices2D[i % end].x, height / 2,
                vertices2D[i % end].y);

            uvs[writeIndex] = center +
                new Vector2(Mathf.Cos(currentRadians), Mathf.Sin(currentRadians)) / 2;

            ++writeIndex;
        }

        return writeIndex;
    }

    /// <summary>
    ///   Updates things and marks as not dirty
    /// </summary>
    private void Update()
    {
        Dirty = false;
        InitializeMesh();
        ApplyAllMaterialParameters();
    }

    // Vector2 GetMovementForCellWall(Vector2 target, Vector2 closestOrganelle);

    /// <summary>
    ///   Cheaper version of contains for absorbing stuff.Calculates a
    ///   circle radius that contains all the points (when it is
    ///   placed at 0,0 local coordinate).
    /// </summary>
    private float CalculateEncompassingCircleRadius()
    {
        if (Dirty)
            Update();

        float distanceSquared = 0;

        foreach (var vertex in vertices2D)
        {
            var currentDistance = vertex.LengthSquared();
            if (currentDistance > distanceSquared)
                distanceSquared = currentDistance;
        }

        return Mathf.Sqrt(distanceSquared);
    }

    private void ApplyAllMaterialParameters()
    {
        ApplyWiggly();
        ApplyMovementWiggly();
        ApplyHealth();
        ApplyTint();
        ApplyTextures();
        ApplyDissolveEffect();
    }

    private void ApplyWiggly()
    {
        if (MaterialToEdit == null)
            return;

        // Don't apply wigglyness too early if this is dirty as getting the circle radius forces membrane position
        // calculation, which we don't want to do twice when initializing a microbe
        if (Dirty)
            return;

        float wigglyNessToApply =
            WigglyNess / (EncompassingCircleRadius * sizeWigglyNessDampeningFactor);

        MaterialToEdit.SetShaderParam("wigglyNess", Mathf.Min(WigglyNess, wigglyNessToApply));
    }

    private void ApplyMovementWiggly()
    {
        if (MaterialToEdit == null)
            return;

        // See comment in ApplyWiggly
        if (Dirty)
            return;

        float wigglyNessToApply =
            MovementWigglyNess / (EncompassingCircleRadius * sizeMovementWigglyNessDampeningFactor);

        MaterialToEdit.SetShaderParam("movementWigglyNess", Mathf.Min(MovementWigglyNess, wigglyNessToApply));
    }

    private void ApplyHealth()
    {
        MaterialToEdit?.SetShaderParam("healthFraction", HealthFraction);
    }

    private void ApplyTint()
    {
        MaterialToEdit?.SetShaderParam("tint", Tint);
    }

    private void ApplyTextures()
    {
        // We must update the texture on already-existing membranes, due to the membrane texture changing
        // for the player microbe.
        if (albedoTexture != null && currentlyLoadedAlbedoTexture == Type.AlbedoTexture)
            return;

        albedoTexture = Type.LoadedAlbedoTexture;

        MaterialToEdit!.SetShaderParam("albedoTexture", albedoTexture);
        MaterialToEdit.SetShaderParam("normalTexture", Type.LoadedNormalTexture);
        MaterialToEdit.SetShaderParam("damagedTexture", Type.LoadedDamagedTexture);
        MaterialToEdit.SetShaderParam("dissolveTexture", noiseTexture);

        currentlyLoadedAlbedoTexture = Type.AlbedoTexture;
    }

    private void ApplyDissolveEffect()
    {
        MaterialToEdit?.SetShaderParam("dissolveValue", DissolveEffectValue);
    }

    private void CopyMeshFromCache(ComputedMembraneData cached)
    {
        // TODO: check if it would be better for us to just keep readonly data in the membrane cache so we could
        // just copy a reference here
        vertices2D.Clear();
        vertices2D.AddRange(cached.Vertices2D);

        // Apply the mesh to us
        Mesh = cached.GeneratedMesh;
        SetSurfaceMaterial(cached.SurfaceIndex, MaterialToEdit);
    }

    /// <summary>
    ///   Creates the actual mesh object. Call InitializeMesh instead of this directly
    /// </summary>
    private void BuildMesh()
    {
        // This is actually a triangle list, but the index buffer is used to build
        // the indices (to emulate a triangle fan)
        var bufferSize = vertices2D.Count + 2;
        var indexSize = vertices2D.Count * 3;

        var arrays = new Array();
        arrays.Resize((int)Mesh.ArrayType.Max);

        // Build vertex, index, and uv lists

        // Index mapping to build all triangles
        var indices = new int[indexSize];
        int currentVertexIndex = 1;

        for (int i = 0; i < indexSize; i += 3)
        {
            indices[i] = 0;
            indices[i + 1] = currentVertexIndex + 1;
            indices[i + 2] = currentVertexIndex;

            ++currentVertexIndex;
        }

        // Write mesh data //
        var vertices = new Vector3[bufferSize];
        var uvs = new Vector2[bufferSize];

        int writeIndex = 0;
        writeIndex = InitializeCorrectMembrane(writeIndex, vertices, uvs);

        if (writeIndex != bufferSize)
            throw new Exception("Membrane buffer write ended up at wrong index");

        // Godot might do this automatically
        // // Set the bounds to get frustum culling and LOD to work correctly.
        // // TODO: make this more accurate by calculating the actual extents
        // m_mesh->_setBounds(Ogre::Aabb(Float3::ZERO, Float3::UNIT_SCALE * 50)
        //     /*, false*/);
        // m_mesh->_setBoundingSphereRadius(50);

        arrays[(int)Mesh.ArrayType.Vertex] = vertices;
        arrays[(int)Mesh.ArrayType.Index] = indices;
        arrays[(int)Mesh.ArrayType.TexUv] = uvs;

        // Create the mesh
        var generatedMesh = new ArrayMesh();

        var surfaceIndex = generatedMesh.GetSurfaceCount();
        generatedMesh.AddSurfaceFromArrays(Mesh.PrimitiveType.Triangles, arrays);

        // Apply the mesh to us
        Mesh = generatedMesh;
        SetSurfaceMaterial(surfaceIndex, MaterialToEdit);

        ProceduralDataCache.Instance.WriteMembraneData(CreateDataForCache(generatedMesh, surfaceIndex));
    }

    private void DrawCorrectMembrane(float cellDimensions, List<Vector2> sourceBuffer, List<Vector2> targetBuffer)
    {
        if (Type.CellWall)
        {
            DrawMembrane(cellDimensions, sourceBuffer, targetBuffer, GetMovementForCellWall);
        }
        else
        {
            DrawMembrane(cellDimensions, sourceBuffer, targetBuffer, GetMovement);
        }
    }

    private void DrawMembrane(float cellDimensions, List<Vector2> sourceBuffer, List<Vector2> targetBuffer,
        Func<Vector2, Vector2, Vector2> movementFunc)
    {
        var tim1 = Time.GetTicksUsec();
        while (targetBuffer.Count < sourceBuffer.Count)
            targetBuffer.Add(new Vector2(0, 0));

        // TODO: check that this is actually needed, and if triggered does the right thing
        while (targetBuffer.Count > sourceBuffer.Count)
            targetBuffer.RemoveAt(targetBuffer.Count - 1);
        var tim2 = Time.GetTicksUsec();
        //// Loops through all the points in the membrane and relocates them as necessary.
        //for (int i = 0, end = sourceBuffer.Count; i < end; ++i)
        //{
        //    var tim11 = Time.GetTicksUsec();
        //    var closestOrganelle = FindClosestOrganelles(sourceBuffer[i]);
        //    var tim12 = Time.GetTicksUsec();
        //    if (closestOrganelle ==
        //        new Vector2(INVALID_FOUND_ORGANELLE, INVALID_FOUND_ORGANELLE))
        //    {
        //        targetBuffer[i] = (sourceBuffer[(end + i - 1) % end] + sourceBuffer[(i + 1) % end]) / 2;
        //    }
        //    else
        //    {
        //        var movementDirection = movementFunc(sourceBuffer[i], closestOrganelle);

        //        targetBuffer[i] = new Vector2(sourceBuffer[i].x - movementDirection.x,
        //            sourceBuffer[i].y - movementDirection.y);
        //    }
        //    var tim13 = Time.GetTicksUsec();
        //    //GD.Print("first " + (tim12 - tim11) + ", second " + (tim13 - tim12));
        //}

        membraneResolution *= 10; //10 works too

        // Loops through all the points in outer square and moves them close to organelles.
        for (int i = 0, end = sourceBuffer.Count; i < end; ++i)
        {
            //var tim11 = Time.GetTicksUsec();
            var closestOrganelle = FindClosestOrganelles(sourceBuffer[i]);
            //var tim12 = Time.GetTicksUsec();
            if (closestOrganelle ==
                new Vector2(INVALID_FOUND_ORGANELLE, INVALID_FOUND_ORGANELLE))
            {
                targetBuffer[i] = (sourceBuffer[(end + i - 1) % end] + sourceBuffer[(i + 1) % end]) / 2;
                //TODO: idk about this, think about it
            }
            else
            {
                var difference = sourceBuffer[i] - closestOrganelle;
                var direction = difference.Normalized();
                float distance = difference.Length();
                var movement = direction * (distance - 2.0f); // 3.0f for eukaryotes because big organelles

                targetBuffer[i] = new Vector2(sourceBuffer[i].x - movement.x,
                    sourceBuffer[i].y - movement.y);
            }

            //var tim13 = Time.GetTicksUsec();
            //GD.Print("first " + (tim12 - tim11) + ", second " + (tim13 - tim12));
        }

        var tim3 = Time.GetTicksUsec();

        float circumference = 0.0f;
        for (int i = 0; i < targetBuffer.Count - 1; ++i)
        {
            circumference += (targetBuffer[i + 1] - targetBuffer[i]).Length();
        }

        // Allows for the addition points in the membrane.
        int pointsAdded = 0;
        for (int i = 0; i < targetBuffer.Count - 1; ++i)
        {
            bool isLast = i == targetBuffer.Count - 1;
            var currentPoint = targetBuffer[i];
            var nextPoint = isLast ? targetBuffer[0] : targetBuffer[i + 1];

            // Check to see if the gap between two points in the membrane is too big.
            if ((nextPoint - currentPoint).Length() >
                circumference / membraneResolution)
            {
                // Add an element after the ith term that is the average of the
                // i and i+1 term.
                //var tempPoint = (targetBuffer[i + 1] + targetBuffer[i]) / 2;

                var direction = (nextPoint - currentPoint).Normalized();
                var newPoint = direction * (circumference / membraneResolution) + targetBuffer[i];
                //var middlePoint = (nextPoint + currentPoint) / 2.0f;

                targetBuffer.Insert(i + 1, newPoint); //TODO: try middle point instead
                pointsAdded++;
                //++i;
            }
        }

        // Allows for the deletion points in the membrane.
        int pointsRemoved = 0;
        for (int i = targetBuffer.Count - 1; i > 2; i--)
        {
            // Check to see if the gap between two points in the membrane is too small.
            if ((targetBuffer[i] - targetBuffer[0]).Length() <
                circumference / membraneResolution)
            {
                // Delete the ith term.
                targetBuffer.RemoveAt(i);
                pointsRemoved++;
            }
            else
            {
                break;
            }
        }

        for (int i = 0; i < targetBuffer.Count - 1; i++)
        {
            // Check to see if the gap between two points in the membrane is too small.
            if ((targetBuffer[i + 1] - targetBuffer[i]).Length() <
                    circumference / membraneResolution)
            {
                // Delete the ith term.
                targetBuffer.RemoveAt(i + 1);
                pointsRemoved++;
                i--;
            }
        }

        // Find the smallest and largest distance between points
        float smallestDistance = float.MaxValue; //squared
        float largestDistance = 0.0f; //squared
        for (int i = 0; i < targetBuffer.Count; ++i)
        {
            float distanceSquared;

            if (i == targetBuffer.Count - 1)
            {
                distanceSquared = (targetBuffer[0] - targetBuffer[i]).LengthSquared();
            }
            else
            {
                distanceSquared = (targetBuffer[i + 1] - targetBuffer[i]).LengthSquared();
            }

            if (distanceSquared < smallestDistance)
            {
                smallestDistance = distanceSquared;
            }

            if (distanceSquared > largestDistance)
            {
                largestDistance = distanceSquared;
            }
        }

        ////Find center
        //var center = new Vector2(0.0f, 0.0f);
        //foreach (var point in targetBuffer)
        //{
        //    center += point;
        //}

        ////Find center 2.0
        //var center = new Vector2(0.0f, 0.0f);
        //foreach (var point in OrganellePositions)
        //{
        //    center += point;
        //}

        //center /= OrganellePositions.Count;

        //var center = combinedVector / targetBuffer.Count;

        //If this isn't a cell wall, make it wavier
        if (!Type.CellWall)
        {
            //Find average distance between center and membrane
            float totalLength = 0.0f;
            for (int i = 0; i < targetBuffer.Count; i++)
            {
                totalLength += targetBuffer[i].Length();//DistanceTo(center);
            }

            float averageLength = totalLength / targetBuffer.Count;

            //Move points a bit to make a less smooth membrane
            float distanceTraveled = 0.0f;
            float multiplier = 2.0f * Mathf.Pi * 10.0f / circumference; //targetBuffer.Count;
            var random = new Random();
            for (int i = 0; i < targetBuffer.Count; i++)
            {
                var point = targetBuffer[i];
                var nextPoint = targetBuffer[(i + 1) % targetBuffer.Count];
                var direction = nextPoint - point;
                distanceTraveled += direction.Length();
                float mmm = Mathf.Sin(multiplier * distanceTraveled) * Mathf.Sqrt(averageLength) * 0.03f;// point.DistanceTo(new Vector2(0.0f, 0.0f)) * 0.1f;
                //Turn 90 degrees
                var newDirection = new Vector2(-direction.y, direction.x);
                point += newDirection.Normalized() * mmm; //new Vector2(newX, newY);
                targetBuffer[i] = point;
            }
        }

        ////Move points a bit to make a less smooth membrane
        ////TODO: try moving away/towards normal instead of center
        //float multiplier = Mathf.Pi / 4;
        //var random = new Random();
        //for (int i = 0; i < targetBuffer.Count; i++)
        //{
        //    var point = targetBuffer[i];
        //    float mmm = Mathf.Sin(i * multiplier) * averageLength * 0.03f;// point.DistanceTo(new Vector2(0.0f, 0.0f)) * 0.1f;
        //    //var newX = point.x + (random.NextFloat() - 0.5f) * multiplier;
        //    //var newY = point.y + (random.NextFloat() - 0.5f) * multiplier;
        //    point += (point - center).Normalized() * mmm; //new Vector2(newX, newY);
        //    targetBuffer[i] = point;
        //}

        ////Randomly move points a bit
        //float multiplier = 0.1f;
        //var random = new Random();
        //for (int i = 0; i < targetBuffer.Count; i++)
        //{
        //    var point = targetBuffer[i];
        //    var newX = point.x + (random.NextFloat() - 0.5f) * multiplier;
        //    var newY = point.y + (random.NextFloat() - 0.5f) * multiplier;
        //    point += point.Normalized() * (random.NextFloat() - 0.5f) * multiplier; //new Vector2(newX, newY);
        //    targetBuffer[i] = point;
        //}

        GD.Print("Points added " + pointsAdded + ", points removed " + pointsRemoved);
        GD.Print("Circumference " + circumference + ", smallest distance " + Mathf.Sqrt(smallestDistance) + ", largest distance " + Mathf.Sqrt(largestDistance));

        //// Allows for the addition and deletion of points in the membrane.
        //for (int i = 0; i < targetBuffer.Count - 1; ++i)
        //{
        //    // Check to see if the gap between two points in the membrane is too big.
        //    if ((targetBuffer[i] - targetBuffer[(i + 1) % targetBuffer.Count]).Length() >
        //        cellDimensions / membraneResolution)
        //    {
        //        // Add an element after the ith term that is the average of the
        //        // i and i+1 term.
        //        var tempPoint = (targetBuffer[(i + 1) % targetBuffer.Count] + targetBuffer[i]) / 2;

        //        targetBuffer.Insert(i + 1, tempPoint);
        //        ++i;
        //    }

        //    // Check to see if the gap between two points in the membrane is too small.
        //    if ((targetBuffer[(i + 1) % targetBuffer.Count] -
        //            targetBuffer[(i + targetBuffer.Count - 1) % targetBuffer.Count]).Length() <
        //        cellDimensions / membraneResolution)
        //    {
        //        // Delete the ith term.
        //        targetBuffer.RemoveAt(i);
        //    }
        //}

        var tim4 = Time.GetTicksUsec();
        //GD.Print("first " + (tim2 - tim1) + ", second " + (tim3 - tim2) + ", third " + (tim4 - tim3));
        GD.Print("Moving points took " + (tim3 - tim2) + ", DrawMembrane took " + (tim4 - tim1) + ", targetBuffer " + targetBuffer.Count);
        membraneResolution /= 10;
    }

    private ComputedMembraneData CreateDataForCache(ArrayMesh mesh, int surfaceIndex)
    {
        // Need to copy our data here when caching it as if we get new organelles and change we would pollute the
        // cache entry
        return new ComputedMembraneData(OrganellePositions, Type, new List<Vector2>(vertices2D), mesh, surfaceIndex);
    }
}
